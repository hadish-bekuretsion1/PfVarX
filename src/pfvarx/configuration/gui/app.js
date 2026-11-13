"use strict";
const byId = id => document.getElementById(id);
const editor = byId("editor");
const saveButton = byId("save");
const status = byId("status");
let state, values, timer, generation = 0;

async function request(path, payload) {
    const response = await fetch(path, payload === undefined ? {} : {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify(payload)
    });
    const result = await response.json();
    if (!response.ok) throw new Error(result.error);
    return result;
}

function resolveSchema(spec) {
    const base = spec.anyOf?.find(item => item.type !== "null") || spec;
    return base.$ref ? state.schema.$defs[base.$ref.split("/").pop()] : base;
}

const filePlaceholders = {
    "config.assay.targets_bed": "inputs/targets/",
    "config.assay.primers_tsv": "inputs/targets/",
    "config.samples_tsv": "inputs/metadata/",
    "config.reference.fasta": "inputs/reference/",
    "config.reference.core_mask_bed": "inputs/reference/",
    "config.reference.gff": "inputs/reference/",
    "config.analyses.resistance_catalog": "inputs/metadata/",
};

function draw(parent, object, schema, path = "") {
    for (const [key, original] of Object.entries(schema.properties || {})) {
        const id = path ? path + "." + key : key;
        const spec = resolveSchema(original);
        if (!(key in object) && "default" in original) object[key] = structuredClone(original.default);
        if (spec.properties) {
            const group = document.createElement("fieldset");
            const title = document.createElement("legend");
            title.textContent = key.replaceAll("_", " ");
            group.append(title);
            group.dataset.field = id;
            draw(group, object[key] ||= {}, spec, id);
            parent.append(group);
            continue;
        }
        const row = document.createElement("div");
        row.dataset.field = id;
        const label = document.createElement("label");
        const labels = {
            "config.assay.design": "Assay design", "config.assay.sequencing_type": "Sequencing type",
            "config.assay.amplification_route": "Amplification route",
            "config.assay.targets_bed": "Target regions (BED)", "config.assay.primers_tsv": "Primer definitions (TSV)",
            "config.reference.core_mask_bed": "Core-genome mask (BED)"
        };
        label.textContent = (labels[id] || key.replaceAll("_", " ")) + " ";
        row.append(label);
        parent.append(row);
        const guide = state.field_guide?.[id];
        if (guide || original.description) {
            const help = document.createElement("small");
            help.id = id + ".help";
            help.textContent = guide ? guide.help + " Values or example: " + guide.example : original.description;
            row.append(help);
        }
        if (id === "config.analyses.requested") {
            for (const [da, contract] of Object.entries(state.analyses)) {
                const choice = document.createElement("label");
                const box = document.createElement("input");
                box.type = "checkbox";
                box.value = da;
                box.checked = (object[key] || []).includes(da);
                box.addEventListener("change", () => {
                    object[key] = Array.from(row.querySelectorAll(":checked"), item => item.value);
                    render();
                });
                choice.append(box, da + " · " + contract.name);
                row.append(choice);
            }
            continue;
        }
        if (id === "config.analyses.selection") {
            for (const value of spec.enum) {
                const choice = document.createElement("label");
                const radio = document.createElement("input");
                radio.type = "radio";
                radio.name = id;
                radio.value = value;
                radio.checked = object[key] === value;
                radio.addEventListener("change", () => {
                    object[key] = value;
                    changed(id);
                });
                choice.append(radio, value === "auto" ? "All design-compatible analyses" : "Select analyses");
                row.append(choice);
            }
            continue;
        }
        let input;
        let options = spec.enum || ({
            "config.analyses.spatial_outcome_source": ["DA002"],
            "config.purpose": ["analysis", "benchmarking"]
        })[id];
        if (id === "config.input_type") options = options.filter(value => value === "FASTQ" || value === "BAM");
        if (id === "config.assay.sequencing_type" && object[key] !== "SWGS") options = options.filter(value => value !== "SWGS");
        if (id === "compute.executor") options = options.filter(value => value !== "slurm");
        if (options) {
            input = document.createElement("select");
            const optional = original.anyOf?.some(item => item.type === "null");
            for (const value of optional ? ["", ...options] : options) {
                const names = {
                    local: "Local", hpc: "HPC (single node)", server: "Server", cloud: "Cloud VM",
                    WGS: "Whole-genome sequencing", SELECTIVE_WGS: "Selective Whole Genome Sequencing",
                    TARGETED_NON_UMI: "Targeted without UMI", TARGETED_UMI: "Targeted with UMI",
                    AMPLICON: "Amplicon", DIRECT: "Direct", SWGA: "Selective whole-genome amplification",
                    TBS_UMI: "Targeted barcode with UMI", TBS_NON_UMI: "Targeted barcode without UMI"
                };
                input.append(new Option(names[value] || value || "Unspecified", value));
            }
            input.value = object[key] ?? "";
            const installation = state.installation;
            if (id === "config.input_type" && installation && installation.profile !== "development" && installation.input_type !== "FASTQ") {
                for (const option of input.options) option.disabled = option.value !== installation.input_type;
                const note = document.createElement("small");
                note.textContent = "This workspace was installed for " + installation.input_type +
                    ". This environment supports aligned input.";
                row.append(note);
            }
        } else if (spec.type === "object" || spec.type === "array") {
            input = document.createElement("textarea");
            input.value = JSON.stringify(object[key] ?? {});
        } else {
            input = document.createElement("input");
            input.type = spec.type === "boolean" ? "checkbox" :
                ["integer", "number"].includes(spec.type) ? "number" : "text";
            input.checked = object[key] === true;
            input.value = object[key] ?? "";
            if (input.type === "number") {
                input.step = spec.type === "integer" ? "1" : "any";
                if (spec.minimum !== undefined) input.min = spec.minimum;
                if (spec.maximum !== undefined) input.max = spec.maximum;
            }
            if (spec.pattern && input.type === "text") input.pattern = spec.pattern;
        }
        input.id = id;
        if (guide || original.description) input.setAttribute("aria-describedby", id + ".help");
        label.htmlFor = id;
        label.append(input);
        const requiredMark = document.createElement("span");
        requiredMark.dataset.requiredMark = "";
        requiredMark.textContent = " *";
        requiredMark.setAttribute("aria-hidden", "true");
        requiredMark.hidden = true;
        label.insertBefore(requiredMark, input);
        input.dataset.schemaRequired = String(
            (schema.required || []).includes(key) && !original.anyOf?.some(item => item.type === "null")
        );
        if (id.endsWith(".schema_version")) input.disabled = true;
        if (options && options.length > 4) {
            const filter = document.createElement("input");
            filter.type = "search";
            filter.placeholder = "Filter " + key.replaceAll("_", " ") + " choices";
            filter.setAttribute("aria-label", filter.placeholder);
            filter.addEventListener("input", () => {
                const query = filter.value.trim().toLowerCase();
                for (const option of input.options) {
                    option.hidden = !option.selected && !(option.text + " " + option.value).toLowerCase().includes(query);
                }
            });
            row.append(filter);
        }
        if (filePlaceholders[id]) input.placeholder = filePlaceholders[id];
        input.addEventListener("input", () => {
            input.setCustomValidity("");
            if (input.type === "checkbox") object[key] = input.checked;
            else if (input.type === "number") object[key] = input.value === "" ? null : input.valueAsNumber;
            else if (input.tagName === "TEXTAREA") {
                try { object[key] = JSON.parse(input.value); }
                catch { input.setCustomValidity("Enter valid JSON"); }
            } else {
                object[key] = input.value === "" && original.anyOf ? null : input.value;
            }
            changed(id);
        });
    }
}

function toggle(id, visible) {
    const row = editor.querySelector('[data-field="' + id + '"]');
    if (!row) return;
    row.hidden = !visible;
    for (const input of row.querySelectorAll("input,select,textarea,button")) input.disabled = !visible;
}

function showInputFolders(directories) {
    const contents = {
        "inputs/fastq": "Paired reads (.fastq or .fastq.gz), such as sample_R1.fastq.gz and sample_R2.fastq.gz.",
        "inputs/alignments": "Coordinate-sorted BAM files and their indexes (.bai or .csi).",
        "inputs/reference": "Reference FASTA and matching core-genome BED. Include matching GFF when annotation is needed.",
        "inputs/metadata": "Sample registry TSV with input paths and study metadata (country, latitude, longitude, collection_date and population where applicable). Optional separate study metadata and verified marker catalogue TSV. Join study metadata to the configured sample registry by sample_id before validation.",
        "inputs/targets": "Target regions BED and primer definitions TSV matched to the reference.",
        "inputs/truth": "truth.json, confident_regions.bed, and <MONOCLONAL_SAMPLE_ID>.truth.vcf for each scored monoclonal sample.",
        "inputs/partitions": "partition.json, created by ./pfvarx benchmark-partition."
    };
    const preview = byId("input-folder-preview");
    preview.replaceChildren();
    for (const directory of directories) {
        const item = document.createElement("li");
        const path = document.createElement("code");
        path.textContent = directory + "/";
        item.append(path, " — " + contents[directory]);
        preview.append(item);
    }
}

function changed(id) {
    const config = values.config;
    const assay = config.assay;
    const analyses = config.analyses;
    if (id === "config.input_type" && config.input_type !== "FASTQ" && config.final_output === "FO001") {
        config.final_output = "FO002";
        render();
        return;
    }
    const location = byId("raw-data-directory");
    const fastq = config.input_type === "FASTQ";
    location.value = fastq ? "inputs/fastq/" : "inputs/alignments/";
    byId("raw-data-help").textContent = fastq
        ? "Example pair: inputs/fastq/sample_R1.fastq.gz and inputs/fastq/sample_R2.fastq.gz. Place your files here after saving, or register their existing paths."
        : "Example: inputs/alignments/sample.bam and inputs/alignments/sample.bam.bai. Place the BAM and its index here after saving, or register their existing paths.";
    const outputs = byId("config.final_output");
    for (const option of outputs.options) {
        option.disabled = option.value === "FO001" && config.input_type !== "FASTQ";
    }
    const full = config.final_output === "FO004";
    const targeted = assay.design !== "WGS";
    if (id === "config.assay.design") {
        assay.sequencing_type = null;
        assay.amplification_route = targeted ? "AMPLICON" : "DIRECT";
        assay.targets_bed = targeted ? assay.targets_bed || filePlaceholders["config.assay.targets_bed"] : null;
        assay.primers_tsv = targeted ? assay.primers_tsv || filePlaceholders["config.assay.primers_tsv"] : null;
        if (assay.design !== "TARGETED_UMI") assay.read_structure = null;
        analyses.requested = (analyses.requested || []).filter(da => state.analyses[da].designs.includes(assay.design));
        analyses.scales = Object.fromEntries(Object.entries(analyses.scales || {}).filter(([da]) => analyses.requested.includes(da)));
        render();
        return;
    }
    if (id === "config.final_output" || id === "config.analyses.selection" || id === "config.caller.backend") {
        if (!full) analyses.selection = "explicit";
        if (!full || analyses.selection === "auto") {
            analyses.requested = [];
            analyses.scales = {};
        }
        if (full && analyses.selection === "explicit" && !config.reference.gff) {
            config.reference.gff = filePlaceholders["config.reference.gff"];
        }
        if (analyses.selection === "auto" && config.reference.gff === filePlaceholders["config.reference.gff"]) {
            config.reference.gff = null;
        }
        if (config.caller.backend === "native") {
            config.caller.executable = null;
            config.caller.sha256 = null;
        }
        render();
        return;
    }
    analyses.scales = Object.fromEntries(
        Object.entries(analyses.scales || {}).filter(([da]) => (analyses.requested || []).includes(da))
    );
    if (full && (analyses.requested || []).includes("DA007") && !analyses.spatial_outcome_source) {
        analyses.spatial_outcome_source = "DA002";
        const source = byId("config.analyses.spatial_outcome_source");
        if (source) source.value = "DA002";
    }
    toggle("config.analyses", full);
    toggle("config.analyses.requested", full && analyses.selection !== "auto");
    toggle("config.analyses.scales", full && analyses.selection !== "auto");
    toggle("config.reference.gff", full);
    toggle("config.assay.targets_bed", targeted);
    toggle("config.assay.primers_tsv", targeted);
    toggle("config.assay.read_structure", assay.design === "TARGETED_UMI");
    toggle("config.caller.executable", config.caller.backend === "external");
    toggle("config.caller.sha256", config.caller.backend === "external");
    toggle("compute.slurm_partition", false);
    const route = byId("config.assay.amplification_route");
    if (id === "config.assay.sequencing_type" && assay.sequencing_type === "SELECTIVE_WGS") {
        assay.amplification_route = "SWGA";
        route.value = "SWGA";
    }
    for (const option of route.options) option.disabled = targeted !== (option.value === "AMPLICON") || (assay.sequencing_type === "SELECTIVE_WGS" && option.value !== "SWGA");
    const sequencing = byId("config.assay.sequencing_type");
    const allowed = {
        WGS: ["WGS", "SELECTIVE_WGS"], TARGETED_NON_UMI: ["AMPLICON", "TBS_NON_UMI"],
        TARGETED_UMI: ["AMPLICON", "TBS_UMI"]
    }[assay.design];
    for (const option of sequencing.options) option.disabled = option.value !== "" && option.value !== assay.sequencing_type && !allowed.includes(option.value);
    for (const box of editor.querySelectorAll('[data-field="config.analyses.requested"] input')) {
        box.disabled = !full || analyses.selection === "auto" || !state.analyses[box.value].designs.includes(assay.design);
    }
    const conditionalRequirements = {
        "config.assay.targets_bed": targeted,
        "config.assay.primers_tsv": targeted,
        "config.assay.read_structure": assay.design === "TARGETED_UMI",
        "config.reference.gff": full && analyses.selection !== "auto",
        "config.analyses.spatial_outcome_source": full && (analyses.requested || []).includes("DA007"),
        "config.caller.executable": config.caller.backend === "external",
        "config.caller.sha256": config.caller.backend === "external"
    };
    for (const input of editor.querySelectorAll("[data-schema-required]")) {
        input.required = !input.disabled && (conditionalRequirements[input.id] ?? (input.dataset.schemaRequired === "true"));
        input.parentElement.querySelector("[data-required-mark]").hidden = !input.required;
    }
    const summary = structuredClone(values);
    if (full && analyses.selection === "explicit") {
        const requested = analyses.requested || [];
        const prerequisites = requested.flatMap(da => state.analyses[da]?.prerequisites || [])
            .map(da => da === "SOURCE" ? analyses.spatial_outcome_source : da);
        summary.internal_prerequisites = [...new Set(prerequisites.filter(da => da && !requested.includes(da)))];
    }
    byId("summary").textContent = JSON.stringify(summary, null, 2);
    saveButton.disabled = false;
    const folders = ["inputs/reference", "inputs/metadata", fastq ? "inputs/fastq" : "inputs/alignments"];
    if (targeted) folders.push("inputs/targets");
    if (config.purpose === "benchmarking") folders.push("inputs/truth", "inputs/partitions");
    showInputFolders(folders);
    const current = ++generation;
    clearTimeout(timer);
    timer = setTimeout(async () => {
        try {
            const invalid = editor.querySelector(":is(input,select,textarea):invalid");
            if (invalid) throw new Error(invalid.id + ": " + invalid.validationMessage);
            const result = await request("check", {values});
            if (current !== generation) return;
            showInputFolders(result.directories);
            status.textContent = "Ready to save.";
        } catch (error) {
            if (current === generation) status.textContent = error.message;
        }
    }, 200);
}

function render() {
    const rawData = byId("raw-data-location");
    editor.replaceChildren();
    draw(editor, values, state.schema);
    byId("config.input_type").closest("[data-field]").after(rawData);
    rawData.hidden = false;
    changed();
}

async function load() {
    try {
        state = await request("state?purpose=" + encodeURIComponent(byId("purpose").value));
        values = state.values;
        if (!state.input_folders_created) {
            const defaults = [
                [values.config.assay, "targets_bed", "inputs/targets/targets.bed"],
                [values.config.assay, "primers_tsv", "inputs/targets/primers.tsv"],
                [values.config, "samples_tsv", "inputs/metadata/samples.tsv"],
                [values.config.reference, "fasta", "inputs/reference/Pf3D7.fasta"],
                [values.config.reference, "core_mask_bed", "inputs/reference/core_genome.bed"]
            ];
            for (const [group, key, template] of defaults) {
                if (group[key] === template) group[key] = filePlaceholders[group === values.config ? "config." + key : group === values.config.assay ? "config.assay." + key : "config.reference." + key];
            }
        }
        render();
    } catch (error) {
        status.textContent = error.message;
    }
}

saveButton.addEventListener("click", async () => {
    ++generation;
    clearTimeout(timer);
    for (const id of Object.keys(filePlaceholders)) {
        const input = byId(id);
        if (input?.value.trim().endsWith("/") && (input.value !== filePlaceholders[id] || id === "config.analyses.resistance_catalog")) input.setCustomValidity("Add the filename after the folder path.");
    }
    const release = byId("config.reference.release_label");
    if (release && ["", "SET_REFERENCE_RELEASE"].includes(release.value.trim())) {
        release.setCustomValidity("Enter the actual reference release label.");
    }
    const invalid = editor.querySelector(":is(input,select,textarea):invalid");
    if (invalid) {
        status.textContent = "Make sure all required fields are completed.";
        invalid.focus();
        invalid.reportValidity();
        return;
    }
    ++generation;
    clearTimeout(timer);
    saveButton.disabled = true;
    try {
        const result = await request("save", {values, revision: state.revision, replace: byId("replace").checked});
        state.revision = result.revision;
        status.textContent = "Input directories are successfully created";
        byId("result").textContent = "Input directories are successfully created\n\nInput folders:\n" + result.directories.map(path => path + "/").join("\n") +
            "\n\nSaved: " + (result.saved.join(", ") || "already current") +
            (result.backups.length ? "\nBackups: " + result.backups.join(", ") : "") +
            "\n\nPlace and register inputs, then run ./pfvarx validate.";
    } catch (error) {
        status.textContent = error.message;
        const field = error.message.match(/config\/(config|quality|compute)\.yaml: ([\w.]+):/);
        if (field) {
            const input = byId(field[1] + "." + field[2]);
            if (input && !input.disabled) {
                input.focus();
                input.setCustomValidity(error.message);
                input.reportValidity();
            }
        }
    } finally {
        saveButton.disabled = false;
    }
});
editor.addEventListener("submit", event => event.preventDefault());
byId("purpose").addEventListener("change", load);
load();

for (const view of ["configuration", "inputs", "validation", "preparation", "preflight"]) {
    byId(view + "-tab").addEventListener("click", () => {
        for (const panel of ["configuration", "inputs", "validation", "preparation", "preflight"]) {
            byId(panel + "-panel").hidden = panel !== view;
            byId(panel + "-tab").setAttribute("aria-pressed", String(panel === view));
        }
        if (view === "inputs") loadFileTabs();
    });
}

for (const [panel, endpoint, title] of [
    ["validation", "validate-inputs", "Input Files Validation"],
    ["preparation", "prepare-inputs", "Reference Preparation"],
    ["preflight", "preflight", "Workflow Readiness Checks"]
]) {
    byId("run-" + panel).addEventListener("click", async () => {
        const button = byId("run-" + panel);
        const message = byId(panel + "-status");
        const output = byId(panel + "-results");
        button.disabled = true;
        output.textContent = "";
        message.textContent = "Running preflight checks…";
        try {
            const result = await request(endpoint, {});
            output.textContent = result.stages.map(stage =>
                "./pfvarx " + stage.command + " (exit " + stage.exit_code + ")\n" +
                [stage.stdout, stage.stderr].filter(Boolean).join("\n")
            ).join("\n");
            message.textContent = result.passed ? title + " passed." :
                title + " stopped at " + result.stages.at(-1).command + ".";
        } catch (error) {
            message.textContent = error.message;
        } finally {
            button.disabled = false;
        }
    });
}

let transferRole;
async function loadFileTabs() {
    const tabs = byId("file-tabs");
    const destination = byId("file-destination");
    tabs.replaceChildren();
    transferRole = undefined;
    try {
        const {folders, extensions, annotation} = await request("input-folders");
        const labels = {
            fastq: [["FASTQ reads", "fastq"]],
            alignments: [["Alignments and indexes", "alignments"]],
            reference: [["Reference FASTA", "reference"], ["Core BED", "reference"], ["GFF annotation", "reference"]],
            metadata: [["Sample registry", "metadata"], ["Study metadata and catalogue", "metadata"]],
            targets: [["Target BED", "targets"], ["Primer TSV", "targets"]],
            truth: [["Truth files", "truth"]]
        };
        const rawFolder = folders.find(folder => ["inputs/fastq", "inputs/alignments"].includes(folder));
        const orderedFolders = rawFolder ? [rawFolder, ...folders.filter(folder => folder !== rawFolder)] : folders;
        for (const folder of orderedFolders) {
            for (const [label, role] of labels[folder.split("/").at(-1)] || []) {
                if (label === "GFF annotation" && !annotation) continue;
                const button = document.createElement("button");
                button.type = "button";
                button.textContent = label;
                button.setAttribute("aria-pressed", "false");
                button.addEventListener("click", () => {
                    transferRole = role;
                    byId("raw-files").value = "";
                    for (const tab of tabs.children) tab.setAttribute("aria-pressed", String(tab === button));
                    const allowed = extensions[role];
                    byId("raw-files").accept = allowed.join(",");
                    destination.textContent = label + " → " + folder + "/. Accepted extensions: " + allowed.join(", ");
                });
                tabs.append(button);
            }
        }
        if (tabs.firstElementChild) tabs.firstElementChild.click();
        else destination.textContent = "Save the configuration and create input directories first.";
    } catch (error) {
        destination.textContent = error.message;
    }
}

async function transfer(action) {
    const controls = byId("inputs-panel").querySelectorAll("button");
    const message = byId("transfer-status");
    controls.forEach(button => { button.disabled = true; });
    message.textContent = "Transferring files…";
    const saved = [];
    try {
        if (!transferRole) throw new Error("Select an input file category first.");
        await action(path => {
            saved.push(path);
            message.textContent = "Saved:\n" + saved.join("\n");
        });
    } catch (error) {
        message.textContent = [...saved.map(path => "Saved: " + path), error.message].join("\n");
    } finally {
        controls.forEach(button => { button.disabled = false; });
    }
}

byId("upload-raw").addEventListener("click", () => transfer(async report => {
    const files = byId("raw-files").files;
    if (!files.length) throw new Error("Choose at least one raw data file.");
    for (const file of files) {
        const response = await fetch("upload-raw?filename=" + encodeURIComponent(file.name) + "&role=" + encodeURIComponent(transferRole), {
            method: "POST", headers: {"Content-Type": "application/octet-stream"}, body: file
        });
        const result = await response.json();
        if (!response.ok) throw new Error(result.error);
        report(result.path);
    }
    byId("raw-files").value = "";
}));

byId("download-raw").addEventListener("submit", event => {
    event.preventDefault();
    const values = Object.fromEntries(new FormData(event.currentTarget));
    transfer(async report => {
        const result = await request("download-raw", {...values, role: transferRole});
        report(result.path + (result.status === "VERIFIED" ? " (checksum verified)" : " (no trusted checksum supplied)") + "\nSHA-256: " + result.sha256);
    });
});
