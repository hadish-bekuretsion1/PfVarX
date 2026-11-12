"use strict";
let D = null;
const $ = (selector, root = document) => root.querySelector(selector);
const el = (tag, attrs = {}, ...children) => { const node = document.createElement(tag); for (const [key, value] of Object.entries(attrs)) { if (key === "class") node.className = value; else if (key.startsWith("on")) node.addEventListener(key.slice(2), value); else if (value !== null && value !== undefined) node.setAttribute(key, value); } for (const child of children.flat(Infinity)) node.append(child instanceof Node ? child : document.createTextNode(child == null ? "" : String(child))); return node; };
const esc = value => String(value ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const number = value => value === null || value === undefined || value === "" || typeof value === "boolean" || (typeof value === "string" && !value.trim()) ? null : (Number.isFinite(Number(value)) ? Number(value) : null);
const fmt = (value, digits = 3) => value === null || value === undefined || value === "" ? "—" : (typeof value === "number" ? (Number.isFinite(value) ? (Number.isInteger(value) ? String(value) : value.toFixed(digits)) : "—") : String(value));
const percent = value => number(value) === null ? "—" : (number(value) * 100).toFixed(1) + "%";
const badge = state => el("span", {class:"badge " + String(state || "UNKNOWN")}, state || "UNKNOWN");
const unique = values => [...new Set(values.filter(v => v !== null && v !== undefined && v !== ""))];
const plain = value => typeof value === "object" && value !== null ? JSON.stringify(value) : (value ?? "");
const note = text => el("div", {class:"note"}, text);
const empty = text => el("div", {class:"empty"}, text || "No values are available in the supplied report data.");
const STATE = {selectedSample:null, redact:false, views:{}, filters:{}, section:"overview", rendered:new Set()};
function notify(message, failed = false) { const box = $("#exportStatus"); box.textContent = message; box.classList.toggle("error", failed); box.style.display = "block"; clearTimeout(notify.timer); notify.timer = setTimeout(() => box.style.display = "none", 6500); }
function showTip(event, html) { const box = $("#tip"); box.innerHTML = html; box.style.display = "block"; const bounds = event.currentTarget?.getBoundingClientRect?.() || {left:20,top:20}; const x = event.clientX || bounds.left + 10, y = event.clientY || bounds.top + 10; box.style.left = Math.max(6, Math.min(x + 12, window.innerWidth - box.offsetWidth - 8)) + "px"; box.style.top = Math.max(6, Math.min(y + 12, window.innerHeight - box.offsetHeight - 8)) + "px"; }
function hideTip() { $("#tip").style.display = "none"; }
function interactive(node, html, action) { node.setAttribute("tabindex", "0"); node.setAttribute("role", action ? "button" : "img"); node.setAttribute("aria-label", html.replace(/<[^>]*>/g, " ")); node.addEventListener("mousemove", e => showTip(e, html)); node.addEventListener("mouseleave", hideTip); node.addEventListener("focus", e => showTip(e, html)); node.addEventListener("blur", hideTip); if (action) { node.addEventListener("click", action); node.addEventListener("keydown", e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); action(); } }); } return node; }
function validateData(input) {
  const errors = [], warnings = [];
  if (!input || typeof input !== "object" || Array.isArray(input)) return {errors:["The report must be a JSON object."], warnings, data:null};
  if (input.schema !== "pfvarx.dashboard/2") errors.push("Expected schema pfvarx.dashboard/2.");
  for (const key of ["project_id", "run_id"]) if (typeof input[key] !== "string" || !input[key].trim()) errors.push("Missing required string: " + key);
  if (!Array.isArray(input.samples)) errors.push("Required samples must be an array.");
  else if (input.samples.some(s => !s || typeof s !== "object" || typeof s.sample_id !== "string")) errors.push("Every sample needs a sample_id string.");
  if (errors.length) return {errors,warnings,data:null};
  const data = {...input};
  for (const key of ["variants","qc_matrix","eligibility","diagnostics","log_content","requested_analyses","internal_analyses"]) { if (data[key] !== undefined && !Array.isArray(data[key])) warnings.push(key + " is invalid; this optional section is unavailable."); data[key] = Array.isArray(data[key]) ? data[key].filter(v => v !== null && v !== undefined) : []; }
  for (const key of ["downstream","resources","software","inputs","run_status","reference","variant_summary"]) { if (!data[key] || typeof data[key] !== "object" || Array.isArray(data[key])) { if (data[key] !== undefined) warnings.push(key + " is invalid; optional details are unavailable."); data[key] = {}; } }
  const record = v => v && typeof v === "object" && !Array.isArray(v);
  for (const key of ["qc_matrix","eligibility","diagnostics","log_content"]) data[key] = data[key].filter(record);
  for (const [id, entry] of Object.entries(data.downstream)) {
    if (!record(entry)) { warnings.push(id + " has invalid optional analysis data; no results substituted."); delete data.downstream[id]; continue; }
    for (const field of ["calls","prevalence","fws","diversity","clonality","pairs","pca","fst","distance_matrix","summary","grid"]) {
      if (entry[field] !== undefined && !Array.isArray(entry[field])) { warnings.push(id + "." + field + " must be an array; unavailable."); entry[field] = []; }
      if (Array.isArray(entry[field])) entry[field] = entry[field].filter(record);
    }
    if (entry.manifest !== undefined && entry.manifest !== null && !record(entry.manifest)) { warnings.push(id + " manifest is invalid; unavailable."); entry.manifest = null; }
  }
  if (data.annotation !== undefined && data.annotation !== null && !record(data.annotation)) { warnings.push("Optional annotation is invalid; unavailable."); data.annotation = null; }
  data.samples = data.samples.map(s => ({...s, metrics:s.metrics && typeof s.metrics === "object" ? s.metrics : {}, permitted:s.permitted || {}, uncalibrated_gates:Array.isArray(s.uncalibrated_gates) ? s.uncalibrated_gates : []}));
  data.variants = data.variants.filter(v => v && typeof v === "object");
  data.resources.measured = Array.isArray(data.resources.measured) ? data.resources.measured.filter(record) : [];
  data.software.parameters = data.software.parameters || {};
  data.software.tools_installed = data.software.tools_installed || {};
  data.report_warnings = [...(Array.isArray(data.report_warnings) ? data.report_warnings : []), ...warnings];
  return {errors,warnings,data};
}
function restoreViews() { try { STATE.views = JSON.parse(sessionStorage.getItem("pfvarx.views." + D.run_id) || "{}"); } catch (_) { STATE.views = {}; } }
function saveViews() { try { sessionStorage.setItem("pfvarx.views." + D.run_id, JSON.stringify(STATE.views)); } catch (_) {} }
function recordsTable(rows, options = {}) { if (!rows.length) return empty("No records match this view or its filters."); const columns = options.columns || unique(rows.flatMap(row => Object.keys(row))).filter(k => !k.startsWith("_")); return table(columns, rows.map(r => columns.map(c => plain(r[c]))), options); }
function numericRows(rows, field) { return rows.filter(r => number(r[field]) !== null); }
function groupBy(rows, key) { const result = new Map(); rows.forEach(r => { const k = String(r[key] ?? "Not recorded"); if (!result.has(k)) result.set(k, []); result.get(k).push(r); }); return result; }
function categoryItems(rows, key) { return [...groupBy(rows, key)].map(([label, list]) => ({label, value:list.length, tip:esc(label) + ": " + list.length + " records"})); }
function visualGroup({id,title,rows,filters = [],views,defaultView,description = "", defaults = {}, extraFilter = null, extraControls = []}) {
  rows = rows.filter(r => r && typeof r === "object");
  const panel = el("div", {class:"panel visual-group", "data-group":id}, el("h3", {}, title));
  if (description) panel.append(note(description));
  const viewId = "view-" + id.replace(/[^a-z0-9]/gi, "-");
  const select = el("select", {id:viewId, "aria-label":"Choose visualization: " + title}, Object.entries(views).map(([value,v]) => el("option", {value}, v.label)));
  const selected = STATE.views[id]; select.value = selected && views[selected] ? selected : (defaultView && views[defaultView] ? defaultView : Object.keys(views)[0]);
  const controls = el("div", {class:"controls"}, el("label", {for:viewId}, "Choose visualization", select));
  const active = el("div", {class:"active-filters", "aria-live":"polite"});
  const body = el("div", {class:"viz-body", "data-view":select.value});
  const state = STATE.filters[id] || (STATE.filters[id] = {...defaults}); const fields = [];
  for (const field of filters) { const values = unique(rows.map(r => r[field])); if (!values.length) continue; const filter = el("select", {"aria-label":"Filter " + field + ": " + title}, el("option", {value:""}, field + ": all"), values.sort((a,b) => String(a).localeCompare(String(b))).map(v => el("option", {value:String(v)}, v))); filter.value = state[field] || ""; if (!filter.value) delete state[field]; filter.onchange = () => { state[field] = filter.value; draw(); }; fields.push([field,filter]); controls.append(filter); }
  controls.append(...extraControls);
  controls.append(el("button", {class:"btn", onclick:() => { fields.forEach(([key,filter]) => { filter.value = ""; delete state[key]; }); extraControls.forEach(c => { if (c._reset) c._reset(); }); draw(); }}, "Reset filters"));
  const current = () => rows.filter(r => fields.every(([field]) => !state[field] || String(r[field]) === state[field])).filter(r => !extraFilter || extraFilter(r));
  const exportCurrent = format => { try { const nested = body.querySelector(".data-table"); if (nested?._export) nested._export(format,id,state); else exportRecords(id,current(),format,state); } catch (error) { notify("Export failed: " + error.message,true); } };
  controls.append(el("button", {class:"btn", onclick:() => exportCurrent("csv")}, "CSV"), el("button", {class:"btn", onclick:() => exportCurrent("xlsx")}, "XLSX"));
  function draw() { const filtered = current(); const labels = fields.filter(([key]) => state[key]).map(([key]) => key + "=" + state[key]); active.textContent = "Showing " + filtered.length + " of " + rows.length + " source records" + (labels.length ? " · Active filters: " + labels.join("; ") : " · No categorical filters"); body.dataset.view = select.value; body.replaceChildren(); if (!filtered.length) body.append(empty("No matching records. Reset filters or choose another dataset.")); else { try { body.append(views[select.value].render(filtered)); } catch (error) { body.append(empty("This view is unavailable: " + error.message)); } } panel.dataset.filteredCount = String(filtered.length); }
  select.onchange = () => { STATE.views[id] = select.value; saveViews(); draw(); }; panel.append(controls,active,body); panel._refresh = draw; draw(); return panel;
}
function loadReport(input) { const check = validateData(input); if (check.errors.length) { $("#main").replaceChildren(el("div", {class:"panel schema-error", role:"alert"}, el("h2", {}, "Report data could not be loaded"), ...check.errors.map(e => el("p", {}, e)), note("Choose a valid dashboard_data.json using Open report JSON. No data was substituted."))); $("#nav").replaceChildren(); $("#dataNotice").textContent = "Invalid report data — no scientific results rendered."; return false; } D = check.data; STATE.selectedSample = null; STATE.filters = {}; STATE.rendered = new Set(); restoreViews(); build(); return true; }
window.PfVarX = {validateData, loadReport, get data() { return D; }, get state() { return STATE; }};
