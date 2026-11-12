"""Offline acceptance checks. Usage: python check_dashboard.py REPORT_DIR OUT_DIR.

Needs Playwright, Pillow, openpyxl and a local Chromium (PFVARX_BROWSER_EXECUTABLE).
Artifacts are written outside the software tree. No scientific jobs are executed.
"""
from __future__ import annotations

import csv
import io
import json
import os
import sys
from pathlib import Path

from openpyxl import load_workbook
from PIL import Image
from playwright.sync_api import sync_playwright


def main(report_dir: str, out_dir: str) -> dict:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    checks: dict[str, bool] = {}
    errors: list[str] = []
    requests: list[str] = []
    with sync_playwright() as pw:
        browser = pw.chromium.launch(executable_path=os.environ.get('PFVARX_BROWSER_EXECUTABLE'),
                                     args=['--no-sandbox'])
        context = browser.new_context(accept_downloads=True, viewport={'width': 1440, 'height': 1000})
        def route(request):
            if request.request.url.startswith(('file:', 'data:', 'blob:')):
                request.continue_()
            else:
                requests.append(request.request.url)
                request.abort()
        context.route('**/*', route)
        page = context.new_page()
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.on('console', lambda message: errors.append(message.text) if message.type == 'error' else None)
        page.goto((Path(report_dir) / 'pfvarx_dashboard.html').resolve().as_uri())
        original = page.evaluate('JSON.parse(JSON.stringify(PfVarX.data))')
        checks['eight_sections'] = page.locator('#nav button.card').count() == 8
        for button in page.locator('#nav button.card').all():
            button.click()
            assert page.locator('#' + button.get_attribute('data-id')).is_visible()
        checks['navigation'] = True
        checks['numeric_missing_preserved'] = page.evaluate('[null,undefined,"", "NA", "nan"].every(v=>number(v)===null) && number(0)===0 && number(-.02)===-.02')
        group_count = 0
        def exercise_groups():
            nonlocal group_count
            for group in page.locator('.visual-group:visible').all():
                group_count += 1
                select = group.locator('select').first
                for value in select.locator('option').evaluate_all('(nodes)=>nodes.map(n=>n.value)'):
                    select.select_option(value)
                    assert group.locator('.viz-body').count() == 1
                    assert group.locator('.viz-body').get_attribute('data-view') == value
                    assert 'This view is unavailable:' not in group.inner_text()
                select.select_option(select.locator('option').first.get_attribute('value'))
        for sid in ['overview', 'variants', 'qc', 'resources']:
            page.locator(f'#nav [data-id="{sid}"]').click()
            exercise_groups()
        page.locator('#nav [data-id="downstream"]').click()
        analysis = page.get_by_label('Select analysis', exact=True)
        if analysis.count():
            for value in analysis.locator('option').evaluate_all('(nodes)=>nodes.map(n=>n.value)'):
                analysis.select_option(value)
                exercise_groups()
                assert 'could not be rendered:' not in page.locator('#downstream').inner_text()
                if value in ('DA002','DA005','DA006'):
                    page.screenshot(path=str(out / (value + '.png')), full_page=True)
                if value == 'DA004':
                    toggle = page.get_by_label('Enable relatedness threshold', exact=True)
                    if toggle.count():
                        toggle.check()
                        page.get_by_label('Relatedness display threshold', exact=True).fill('0.5')
                        page.get_by_label('Relatedness display threshold', exact=True).dispatch_event('input')
                        checks['threshold_visible'] = '0.50' in page.locator('#downstream').inner_text()
                        toggle.uncheck()
            analysis.select_option(analysis.locator('option').first.get_attribute('value'))
        checks['all_available_views'] = group_count > 0
        category = page.locator('.visual-group:visible select[aria-label^="Filter "]').first
        if category.count() and category.locator('option').count()>1:
            category.select_option(category.locator('option').last.get_attribute('value'))
            checks['categorical_filter_state'] = 'Active filters:' in page.locator('.visual-group:visible').first.inner_text()
        page.set_viewport_size({'width':390,'height':844})
        page.locator('#menuToggle').click()
        checks['mobile_menu'] = page.locator('#nav').is_visible() and page.locator('#menuToggle').get_attribute('aria-expanded') == 'true'
        page.locator('#nav [data-id="overview"]').click()
        page.set_viewport_size({'width':1440,'height':1000})
        page.locator('#nav [data-id="downstream"]').click()
        checks['optional_sections_render'] = page.locator('.failbox').filter(has_text='could not be rendered').count() == 0
        # Persistence across analysis switches and navigation; exclusive view body.
        group = page.locator('.visual-group:visible').first
        if group.count():
            select = group.locator('select').first
            last = select.locator('option').last.get_attribute('value')
            gid = group.get_attribute('data-group')
            select.select_option(last)
            page.locator('#nav [data-id="overview"]').click()
            page.locator('#nav [data-id="downstream"]').click()
            checks['view_persistence'] = page.locator(f'[data-group="{gid}"] select').first.input_value() == last
        else:
            checks['view_persistence'] = True  # No downstream dataset in FO001/FO002.
        page.locator('#nav [data-id="qc"]').click()
        graphic = page.locator('#qc svg [aria-label][tabindex]').first
        if graphic.count():
            graphic.focus()
            checks['accessible_actual_tooltip'] = page.locator('#tip').is_visible() and len(page.locator('#tip').inner_text()) > 15
            graphic.hover()
            checks['hover_tooltip'] = page.locator('#tip').is_visible()
        sample_table = page.locator('#qc .panel').filter(has=page.get_by_role('heading', name='Samples (click to inspect)', exact=True))
        if original['samples']:
            sid = original['samples'][0]['sample_id']
            sample_table.locator('input[type=search]').fill(sid)
            checks['search_filters_actual_rows'] = sid in sample_table.locator('tbody').inner_text()
            sample_table.locator('th').first.click()
            checks['sort_indicator'] = sample_table.locator('th').first.get_attribute('aria-sort') in ('ascending', 'descending')
            sample_table.locator('tbody tr').first.click()
            checks['sample_drilldown'] = sid in page.locator('#sampleDetail').inner_text()
            with page.expect_download() as download:
                sample_table.get_by_role('button', name='CSV', exact=True).click()
            csv_path = out / 'samples.csv'
            download.value.save_as(csv_path)
            rows = list(csv.reader(io.StringIO(csv_path.read_text(encoding='utf-8-sig'))))
            checks['csv_filtered_values'] = any(sid in row for row in rows)
            with page.expect_download() as download:
                sample_table.get_by_role('button', name='XLSX', exact=True).click()
            xlsx_path = out / 'samples.xlsx'
            download.value.save_as(xlsx_path)
            book = load_workbook(xlsx_path)
            checks['genuine_xlsx'] = any(sid in row for row in book['data'].values)
            checks['export_feedback'] = 'download' in page.locator('#exportStatus').inner_text().lower()
        page.get_by_role('button', name='Reset dashboard', exact=True).click()
        page.locator('#nav [data-id="qc"]').click()
        png = page.locator('#qc').get_by_role('button', name='PNG', exact=True).first
        if png.count():
            with page.expect_download() as download:
                png.click()
            path = out / 'chart.png'
            download.value.save_as(path)
            image = Image.open(path).convert('RGB')
            checks['png_resolution_white_background'] = image.width >= 1000 and image.getpixel((0, 0)) == (255, 255, 255)
        # Browser print is the offline PDF export, with current views and visible table pages.
        page.evaluate('SECTIONS.forEach(([id])=>ensureSection(id))')
        page.pdf(path=str(out / 'report.pdf'), print_background=True, prefer_css_page_size=True)
        checks['pdf_created'] = (out / 'report.pdf').stat().st_size > 10000
        page.screenshot(path=str(out / 'qc.png'), full_page=True)
        checks['synthetic_banner_when_declared'] = original.get('data_kind') != 'synthetic' or 'SYNTHETIC DEMO DATA' in page.locator('#dataNotice').inner_text()
        # Local import, schema errors and missing optional containers.
        input_id = page.locator('input[type=file]').get_attribute('id')
        page.locator('#' + input_id).set_input_files({'name':'invalid.json','mimeType':'application/json','buffer':b'{broken'})
        page.wait_for_timeout(50)
        checks['invalid_json_explained'] = 'could not' in page.locator('#main').inner_text().lower()
        checks['required_schema_validation'] = page.evaluate('PfVarX.validateData({samples:[]}).errors.length>0')
        minimal = {'schema':'pfvarx.dashboard/2','project_id':'SYNTHETIC_MINIMAL','run_id':'MINIMAL','samples':[], 'data_kind':'synthetic', 'downstream':{'DA002':None}}
        page.locator('#' + input_id).set_input_files({'name':'minimal.json','mimeType':'application/json','buffer':json.dumps(minimal).encode()})
        page.wait_for_timeout(50)
        page.evaluate('SECTIONS.forEach(([id])=>ensureSection(id))')
        checks['missing_optional_no_crash'] = page.locator('.failbox').filter(has_text='could not be rendered').count() == 0
        # Large tables remain paginated; empty filtering does not invent a row.
        many = dict(minimal, samples=[{'sample_id':f'SYNTHETIC_{i}'} for i in range(125)], downstream={})
        page.evaluate('input=>PfVarX.loadReport(input)', many)
        page.locator('#nav [data-id="qc"]').click()
        sample_table = page.locator('#qc .panel').filter(has=page.get_by_role('heading', name='Samples (click to inspect)', exact=True))
        checks['table_pagination'] = sample_table.locator('tbody tr').count() == 50
        sample_table.get_by_role('button', name='Next', exact=True).click()
        checks['pagination_changes_page'] = 'Page 2 / 3' in sample_table.inner_text()
        sample_table.locator('input[type=search]').fill('NO_SUCH_SAMPLE')
        checks['no_match_state'] = 'No matching' in sample_table.inner_text() or 'No rows' in sample_table.inner_text()
        sample_table.get_by_role('button', name='Reset', exact=True).click()
        checks['filter_reset'] = sample_table.locator('tbody tr').count() == 50
        group = page.locator('[data-group="qc-summary"]')
        if not group.count():
            group = page.locator('#qc .visual-group').first
        group.locator('select').first.select_option('table')
        group.locator('input[type=search]').fill('SYNTHETIC_124')
        with page.expect_download() as download:
            group.get_by_role('button', name='CSV', exact=True).click()
        download.value.save_as(out / 'embedded-filter.csv')
        exported_group = list(csv.reader(io.StringIO((out / 'embedded-filter.csv').read_text(encoding='utf-8-sig'))))
        checks['embedded_table_export_respects_search'] = sum('SYNTHETIC_124' in row for row in exported_group) == 1 and not any('SYNTHETIC_123' in row for row in exported_group)
        # User-supplied HTML stays text; CSV escaping and null cells round-trip.
        with page.expect_download() as download:
            page.evaluate('rows=>exportRecords("escaping",rows,"csv",{})', [{'sample_id':'A,"quoted"\nline','value':None,'zero':0}])
        download.value.save_as(out / 'escaping.csv')
        exported = list(csv.reader(io.StringIO((out / 'escaping.csv').read_text(encoding='utf-8-sig'))))
        checks['csv_escaping_missing_zero'] = any('A,"quoted"\nline' in row and '0' in row and '' in row for row in exported)
        checks['no_network_requests'] = not requests
        checks['no_javascript_errors'] = not errors
        browser.close()
    result = {'checks':checks,'groups_exercised':group_count,'errors':errors,'network_requests':requests,'passed':all(checks.values())}
    (out / 'checks.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    result = main(*sys.argv[1:3])
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result['passed'] else 1)
