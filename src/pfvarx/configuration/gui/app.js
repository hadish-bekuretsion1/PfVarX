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

function button(label, action) {
    const element = document.createElement("button");
    element.type = "button";
    element.textContent = label;
    element.addEventListener("click", action);
    return element;
}

async function browse(path, input) {
    const picker = byId("picker");
    const list = byId("entries");
    list.replaceChildren();
    if (!picker.open) picker.showModal();
    try {
        const entries = await request("browse?path=" + encodeURIComponent(path));
        if (path !== ".") {
            list.append(button("Parent folder", () => browse(path.split("/").slice(0, -1).join("/") || ".", input)));
            list.append(button("Use this directory", () => selectPath(path, input)));
        }
        const search = document.createElement("input");
        search.type = "search";
        search.placeholder = "Filter files and folders";
        search.setAttribute("aria-label", "Filter files and folders");
        list.append(search);
        const choices = document.createElement("div");
        list.append(choices);
        search.addEventListener("input", () => {
            for (const item of choices.children) {
                item.hidden = !item.textContent.toLowerCase().includes(search.value.toLowerCase());
            }
        });
        for (const entry of entries) {
            choices.append(button(entry.path + (entry.directory ? "/" : ""), () =>
                entry.directory ? browse(entry.path, input) : selectPath(entry.path, input)));
        }
        if (!entries.length) list.append("No entries in this directory.");
    } catch (error) {
        list.textContent = error.message;
    }
}

function selectPath(path, input) {
    input.value = "inputs/" + path;
    input.dispatchEvent(new Event("input"));
    byId("picker").close();
}

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
                box.checked = da === "DA001" || (object[key] || []).includes(da);
                box.addEventListener("change", () => {
                    object[key] = Array.from(row.querySelectorAll(":checked"), item => item.value);
                    changed(id);
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
        if (id === "compute.executor") options = options.filter(value => value !== "slurm");
        if (options) {
            input = document.createElement("select");
            const optional = original.anyOf?.some(item => item.type === "null");
            for (const value of optional ? ["", ...options] : options) {
                const names = {
                    local: "Local", hpc: "HPC (single node)", server: "Server", cloud: "Cloud VM",
                    WGS: "Whole-genome sequencing", SWGS: "Shallow whole-genome sequencing",
                    TARGETED_NON_UMI: "Targeted without UMI", TARGETED_UMI: "Targeted with UMI",
                    AMPLICON: "Amplicon", DIRECT: "Direct", SWGA: "Selective whole-genome amplification",
                    TBS_UMI: "Targeted barcode with UMI", TBS_NON_UMI: "Targeted barcode without UMI"
                };
                input.append(new Option(names[value] || value || "Unspecified", value));
            }
            input.value = object[key] ?? "";
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
        requiredMark.textContent = " (required)";
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
        if (spec.format === "path" && input.type === "text") {
            row.append(button("Browse", () => browse(".", input)));
        }
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

function changed(id) {
    const config = values.config;
    const assay = config.assay;
    const analyses = config.analyses;
    if (id === "config.input_type" && config.input_type !== "FASTQ" && config.final_output === "FO001") {
        config.final_output = "FO002";
        render();
        return;
    }
    const outputs = byId("config.final_output");
    for (const option of outputs.options) {
        option.disabled = option.value === "FO001" && config.input_type !== "FASTQ";
    }
    const full = config.final_output === "FO004";
    const targeted = assay.design !== "WGS";
    if (id === "config.assay.design") {
        assay.sequencing_type = null;
        assay.amplification_route = targeted ? "AMPLICON" : "DIRECT";
        assay.targets_bed = targeted ? assay.targets_bed || "inputs/targets/targets.bed" : null;
        assay.primers_tsv = targeted ? assay.primers_tsv || "inputs/targets/primers.tsv" : null;
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
        if (config.caller.backend === "native") {
            config.caller.executable = null;
            config.caller.sha256 = null;
        }
        render();
        return;
    }
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
    for (const option of route.options) option.disabled = targeted !== (option.value === "AMPLICON");
    const sequencing = byId("config.assay.sequencing_type");
    const allowed = {
        WGS: ["WGS", "SWGS"], TARGETED_NON_UMI: ["AMPLICON", "TBS_NON_UMI"],
        TARGETED_UMI: ["AMPLICON", "TBS_UMI"]
    }[assay.design];
    for (const option of sequencing.options) option.disabled = option.value !== "" && !allowed.includes(option.value);
    for (const box of editor.querySelectorAll('[data-field="config.analyses.requested"] input')) {
        box.disabled = box.value === "DA001" || !full || analyses.selection === "auto" || !state.analyses[box.value].designs.includes(assay.design);
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
    byId("summary").textContent = JSON.stringify(values, null, 2);
    saveButton.disabled = true;
    const current = ++generation;
    clearTimeout(timer);
    timer = setTimeout(async () => {
        try {
            const invalid = editor.querySelector(":invalid:not(form)");
            if (invalid) throw new Error(invalid.id + ": " + invalid.validationMessage);
            const result = await request("check", {values});
            if (current !== generation) return;
            status.textContent = "Ready to save. Input folders: " + result.directories.join(", ");
            saveButton.disabled = false;
        } catch (error) {
            if (current === generation) status.textContent = error.message;
        }
    }, 200);
}

function render() {
    editor.replaceChildren();
    draw(editor, values, state.schema);
    changed();
}

async function load() {
    try {
        state = await request("state?purpose=" + encodeURIComponent(byId("purpose").value));
        values = state.values;
        render();
    } catch (error) {
        status.textContent = error.message;
    }
}

saveButton.addEventListener("click", async () => {
    saveButton.disabled = true;
    try {
        const result = await request("save", {values, revision: state.revision, replace: byId("replace").checked});
        state.revision = result.revision;
        byId("result").textContent = "Input folders:\n" + result.directories.join("\n") +
            "\n\nSaved: " + (result.saved.join(", ") || "already current") +
            (result.backups.length ? "\nBackups: " + result.backups.join(", ") : "") +
            "\n\nPlace and register inputs, then run ./pfvarx validate.";
    } catch (error) {
        status.textContent = error.message;
    } finally {
        saveButton.disabled = false;
    }
});
editor.addEventListener("submit", event => event.preventDefault());
byId("purpose").addEventListener("change", load);
byId("close").addEventListener("click", () => byId("picker").close());
load();
