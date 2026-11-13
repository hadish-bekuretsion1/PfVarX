"""Browser acceptance test for the optional configuration editor."""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from playwright.sync_api import expect, sync_playwright

source = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="pfvarx040-gui-") as directory:
    ws = Path(directory)
    shutil.copytree(source / "config", ws / "config")
    shutil.copytree(source / "workflow", ws / "workflow")
    shutil.copyfile(source / "pixi.toml", ws / "pixi.toml")
    python = os.environ.get("PFVARX_GUI_PYTHON", "python")
    env = dict(os.environ, PYTHONPATH=str(source / "src"), PYTHONDONTWRITEBYTECODE="1")
    process = subprocess.Popen(
        [str(python), "-m", "pfvarx.cli", "configure", "--gui", "--no-browser"],
        cwd=ws,
        env=env,
        stdout=subprocess.PIPE,
        text=True,
    )
    try:
        assert process.stdout
        url = process.stdout.readline().strip().split()[-1]
        assert url.startswith("http://127.0.0.1:"), url
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(
                executable_path=os.environ.get("PFVARX_BROWSER_EXECUTABLE"), args=["--no-sandbox"]
            )
            page = browser.new_page(viewport={"width": 390, "height": 900})
            errors = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(url)
            source_type = page.locator('[id="config.input_type"]')
            output = page.locator('[id="config.final_output"]')
            assert source_type.locator('option[value="VCF"]').count() == 0
            assert source_type.locator("option").all_text_contents() == ["FASTQ", "BAM"]
            for aligned in ("BAM",):
                source_type.select_option("FASTQ")
                output.select_option("FO001")
                source_type.select_option(aligned)
                expect(output).to_have_value("FO002")
                assert output.locator('option[value="FO001"]').evaluate('(option) => option.disabled')
            source_type.select_option("FASTQ")
            assert page.locator('[id="compute.executor"] option').evaluate_all(
                "(options) => options.map(option => option.value)"
            ) == ["local", "hpc", "server", "cloud"]
            save = page.locator("#save")
            expect(save).to_be_disabled()
            page.locator('[id="config.reference.release_label"]').fill("fixture release")
            page.locator('[id="config.assay.design"]').select_option("TARGETED_UMI")
            expect(page.locator('[id="config.assay.read_structure"]')).to_be_visible()
            expect(page.locator('[id="config.assay.amplification_route"]')).to_have_value("AMPLICON")
            expect(save).to_be_disabled()
            page.locator('[id="config.assay.read_structure"]').fill("8M+T 8M+T")
            page.get_by_role("searchbox", name="Filter sequencing type choices").fill("TBS")
            page.locator('[id="config.assay.sequencing_type"]').select_option("TBS_UMI")
            page.locator('[id="config.final_output"]').select_option("FO004")
            page.locator('[name="config.analyses.selection"][value="auto"]').check()
            expect(page.locator('[data-field="config.analyses.requested"]')).to_be_hidden()
            page.locator('[id="config.reference.gff"]').fill("inputs/reference/fixture.gff")
            expect(save).to_be_enabled()
            save.click()
            expect(page.locator("#result")).to_contain_text("inputs/targets")
            assert (ws / "inputs/targets").is_dir()
            assert not (ws / "inputs/metadata/samples.tsv").exists()
            assert not (ws / "prepared").exists()
            save.click()
            expect(page.locator("#result")).to_contain_text("already current")
            (ws / "inputs/reference/fixture.gff").write_text("fixture")
            page.locator('[data-field="config.reference.gff"]').get_by_role("button", name="Browse").click()
            page.get_by_role("button", name="reference/").click()
            page.get_by_role("searchbox", name="Filter files and folders").fill("fixture")
            page.get_by_role("button", name="reference/fixture.gff").click()
            expect(page.locator('[id="config.reference.gff"]')).to_have_value("inputs/reference/fixture.gff")
            assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
            assert not errors, errors
            print("GUI passed validation, exclusive choices, conditional fields, search, save and idempotent rerun")
    finally:
        process.terminate()
        process.wait(timeout=10)
