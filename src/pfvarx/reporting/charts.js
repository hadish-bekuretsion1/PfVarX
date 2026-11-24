const NS = "http://www.w3.org/2000/svg";
const S = (tag, a = {}) => { const e = document.createElementNS(NS, tag); for (const [k, v] of Object.entries(a)) e.setAttribute(k, v); return e; };
function zoomable(svg, w, h) { let vb = [0, 0, w, h], drag = null;
  svg.addEventListener("wheel", e => { e.preventDefault(); const f = e.deltaY > 0 ? 1.15 : 1 / 1.15; const r = svg.getBoundingClientRect(); const mx = vb[0] + (e.clientX - r.left) / r.width * vb[2], my = vb[1] + (e.clientY - r.top) / r.height * vb[3];
    vb = [mx - (mx - vb[0]) * f, my - (my - vb[1]) * f, vb[2] * f, vb[3] * f]; svg.setAttribute("viewBox", vb.join(" ")); }, { passive: false });
  svg.addEventListener("mousedown", e => drag = [e.clientX, e.clientY, vb[0], vb[1]]); svg.addEventListener("mouseup", () => drag = null); svg.addEventListener("mouseleave", () => drag = null);
  svg.addEventListener("mousemove", e => { if (!drag) return; const r = svg.getBoundingClientRect(); vb[0] = drag[2] - (e.clientX - drag[0]) / r.width * vb[2]; vb[1] = drag[3] - (e.clientY - drag[1]) / r.height * vb[3]; svg.setAttribute("viewBox", vb.join(" ")); });
  return () => { vb = [0, 0, w, h]; svg.setAttribute("viewBox", vb.join(" ")); }; }
function chartPanel(title,svg,name,extra = []) {
  if (!(svg instanceof SVGElement)) return el("div",{class:"chart-panel"},svg);
  svg.setAttribute("font-family","system-ui, sans-serif");
  for (const node of svg.querySelectorAll("[tabindex],circle,path")) {
    node.setAttribute("tabindex","0"); node.setAttribute("aria-describedby","tip");
    node.addEventListener("focus",() => { hideTip(); const bounds=node.getBoundingClientRect(); node.dispatchEvent(new MouseEvent("mousemove",{clientX:bounds.left+4,clientY:bounds.top+4})); if ($("#tip").style.display === "block") node.setAttribute("aria-label",$("#tip").textContent); });
    node.addEventListener("blur",hideTip);
    node.addEventListener("keydown",event => { if(event.key === "Enter" || event.key === " ") {event.preventDefault();node.dispatchEvent(new MouseEvent("click",{bubbles:true}));} });
  }
  const reset = svg._reset;
  return el("div",{class:"chart-panel"},el("div",{class:"controls"},...extra,reset?el("button",{class:"btn",onclick:reset},"Reset zoom"):null,el("button",{class:"btn",onclick:()=>exportPNG(svg,name,title)},"PNG")),el("div",{class:"chart-host"},svg),svg._note?note(svg._note):null);
}
function axisText(svg,text,x,y,attrs = {}) { const label=S("text",{x,y,"font-size":11,fill:"#455d70",...attrs});label.textContent=text;svg.append(label);return label; }
function barChart(items,opts = {}) {
  if (!items.length) return empty("No values are available for this chart.");
  const shown=items.slice(0,60),w=760,h=340,m={l:65,r:16,t:20,b:105}; const values=shown.flatMap(i=>[number(i.value),number(i.lo),number(i.hi)]).filter(v=>v!==null);
  if(!values.length)return empty("All values in this chart are missing; no zero values were substituted.");
  const low=Math.min(0,...values),high=opts.max ?? Math.max(...values,1e-9);const range=high-low||1;const y=value=>h-m.b-(value-low)/range*(h-m.t-m.b);const baseline=y(0),bw=(w-m.l-m.r)/shown.length;
  const svg=S("svg",{viewBox:`0 0 ${w} ${h}`,width:w,role:"img","aria-label":opts.label||"bar chart"});
  for(let k=0;k<=4;k++){const value=low+range*k/4;svg.append(S("line",{x1:m.l,y1:y(value),x2:w-m.r,y2:y(value),stroke:"#e4ebf1"}));axisText(svg,fmt(value,2),m.l-8,y(value)+4,{"text-anchor":"end"});}
  shown.forEach((item,i)=>{const x=m.l+i*bw+bw*.15,value=number(item.value);if(value===null)axisText(svg,"missing",x+2,h-m.b-7,{"font-size":9});else{const rect=S("rect",{x,y:Math.min(y(value),baseline),width:Math.max(1,bw*.7),height:Math.max(1,Math.abs(baseline-y(value))),fill:item.color||"#27749d"});interactive(rect,item.tip||`${esc(item.label)}: ${fmt(value)}${number(item.lo)!==null?"; interval "+fmt(item.lo)+"–"+fmt(item.hi):""}`,opts.onClick?()=>opts.onClick(item):null);svg.append(rect);if(number(item.lo)!==null&&number(item.hi)!==null){const center=x+bw*.35;svg.append(S("line",{x1:center,x2:center,y1:y(number(item.lo)),y2:y(number(item.hi)),stroke:"#17364b","stroke-width":1.5}));}}
    axisText(svg,String(item.label),x+bw*.35,h-m.b+16,{"font-size":10,"text-anchor":"end",transform:`rotate(-45 ${x+bw*.35} ${h-m.b+16})`});});
  axisText(svg,opts.ylab||"Value",13,110,{transform:"rotate(-90 13 110)","text-anchor":"middle"});svg._reset=zoomable(svg,w,h);if(items.length>shown.length)svg._note=`First ${shown.length} of ${items.length} categories shown; the table and exports retain all filtered records.`;return svg;
}
function forestPlot(items,opts = {}) {
  const w=850,labels=300,shown=items.slice(0,40),h=80+shown.length*32;const svg=S("svg",{viewBox:`0 0 ${w} ${h}`,width:w,role:"img","aria-label":"forest plot"});const valid=shown.flatMap(r=>[number(r.value),number(r.lo),number(r.hi)]).filter(v=>v!==null);if(!valid.length)return empty("No supplied point estimates or confidence intervals are available.");const min=opts.min??Math.min(0,...valid),max=opts.max??Math.max(...valid,1e-9),x=v=>labels+(v-min)/(max-min||1)*(w-labels-145);
  for(let k=0;k<=4;k++){const value=min+(max-min)*k/4;svg.append(S("line",{x1:x(value),x2:x(value),y1:15,y2:h-43,stroke:"#e4ebf1"}));axisText(svg,fmt(value,2),x(value),h-24,{"text-anchor":"middle"});}
  shown.forEach((r,i)=>{const y=34+i*32;axisText(svg,r.label,labels-13,y+4,{"text-anchor":"end","font-size":11});if(number(r.value)===null){axisText(svg,"Not available",labels+5,y+4);return;}if(number(r.lo)!==null&&number(r.hi)!==null)svg.append(S("line",{x1:x(number(r.lo)),x2:x(number(r.hi)),y1:y,y2:y,stroke:"#337b9c","stroke-width":2}));const point=S("circle",{cx:x(number(r.value)),cy:y,r:5,fill:"#17628b"});interactive(point,r.tip||`${esc(r.label)}: ${fmt(r.value)}; interval ${fmt(r.lo)}–${fmt(r.hi)}`,opts.onClick?()=>opts.onClick(r):null);svg.append(point);axisText(svg,r.count||"",w-125,y+4,{"font-size":10});});axisText(svg,opts.xlab||"Estimate and supplied interval",labels+(w-labels-145)/2,h-3,{"text-anchor":"middle"});if(items.length>shown.length)svg._note=`First ${shown.length} estimates shown; use filters or the table for the full dataset.`;svg._reset=zoomable(svg,w,h);return svg;
}
function histogram(values,bins=20,label="histogram",opts={}) {
  const valid=values.map(number).filter(v=>v!==null);if(!valid.length)return empty("No observed numerical values are available for this distribution.");let min=opts.min??valid.reduce((a,b)=>Math.min(a,b),Infinity),max=opts.max??valid.reduce((a,b)=>Math.max(a,b),-Infinity);if(min===max){min-=.5;max+=.5;}const counts=Array(bins).fill(0);for(const value of valid){const index=Math.max(0,Math.min(bins-1,Math.floor((value-min)/(max-min)*bins)));counts[index]++;}
  const chart=barChart(counts.map((count,i)=>{const low=min+i*(max-min)/bins,high=min+(i+1)*(max-min)/bins;return {label:fmt(low,2)+"–"+fmt(high,2),value:count,tip:`${esc(label)}: ${fmt(low)}–${fmt(high)}; count ${count}; ${(100*count/valid.length).toFixed(2)}% of ${valid.length} observed values`};}),{label,ylab:"Number of observed values"});return chart;
}
function scatter(points,opts={}) {
  const valid=points.filter(p=>number(p.x)!==null&&number(p.y)!==null);if(!valid.length)return empty("Both measured axes are required for this scatter plot.");const shown=valid.slice(0,600),w=760,h=450,m={l:75,r:170,t:32,b:62};let x0=Infinity,x1=-Infinity,y0=Infinity,y1=-Infinity;for(const p of valid){x0=Math.min(x0,number(p.x));x1=Math.max(x1,number(p.x));y0=Math.min(y0,number(p.y));y1=Math.max(y1,number(p.y));}const sx=v=>m.l+(v-x0)/(x1-x0||1)*(w-m.l-m.r),sy=v=>h-m.b-(v-y0)/(y1-y0||1)*(h-m.t-m.b);const svg=S("svg",{viewBox:`0 0 ${w} ${h}`,width:w,role:"img","aria-label":opts.label||"scatter"});const groups=unique(valid.map(p=>p.group??"Not recorded")),colors=["#24759e","#a65f31","#4c8c69","#7a65a3","#ba7a33","#438b8e"];
  for(let k=0;k<=4;k++){const xx=x0+(x1-x0||1)*k/4,yy=y0+(y1-y0||1)*k/4;svg.append(S("line",{x1:m.l,x2:w-m.r,y1:sy(yy),y2:sy(yy),stroke:"#e4ebf1"}));axisText(svg,fmt(xx,2),sx(xx),h-m.b+20,{"text-anchor":"middle"});axisText(svg,fmt(yy,2),m.l-8,sy(yy)+4,{"text-anchor":"end"});}
  shown.forEach(p=>{const node=S("circle",{cx:sx(number(p.x)),cy:sy(number(p.y)),r:STATE.selectedSample===p.label?7:4.8,fill:colors[groups.indexOf(p.group??"Not recorded")%colors.length],stroke:"white"});interactive(node,p.tip||`${esc(p.label)}; ${esc(opts.xlab||"x")} ${fmt(p.x)}; ${esc(opts.ylab||"y")} ${fmt(p.y)}; group ${esc(p.group??"Not recorded")}`,opts.onClick?()=>opts.onClick(p):null);svg.append(node);});groups.slice(0,20).forEach((group,i)=>{svg.append(S("rect",{x:w-m.r+20,y:20+i*17,width:9,height:9,fill:colors[i%colors.length]}));axisText(svg,group,w-m.r+35,29+i*17,{"font-size":10});});axisText(svg,opts.xlab||"x",(w-m.r+m.l)/2,h-8,{"text-anchor":"middle"});axisText(svg,opts.ylab||"y",15,h/2,{transform:`rotate(-90 15 ${h/2})`,"text-anchor":"middle"});svg._reset=zoomable(svg,w,h);if(shown.length<valid.length)svg._note=`First ${shown.length} of ${valid.length} observed points drawn; exports retain the full filtered data.`;return svg;
}
function heatmap(rowsL,colsL,value,opts={}) {
  if(!rowsL.length||!colsL.length)return empty("No rows or columns are available for this heatmap.");const cols=colsL.slice(0,40),rows=rowsL.slice(0,Math.min(50,Math.max(1,Math.floor(2000/cols.length))));const cw=24,ch=22,ml=180,mt=145,w=ml+cols.length*cw+15,h=mt+rows.length*ch+15;const svg=S("svg",{viewBox:`0 0 ${w} ${h}`,width:Math.min(w,1150),role:"img","aria-label":opts.label||"heatmap"});cols.forEach((col,j)=>axisText(svg,col,ml+j*cw+14,mt-8,{"font-size":10,transform:`rotate(-55 ${ml+j*cw+14} ${mt-8})`}));rows.forEach((row,i)=>{axisText(svg,row,ml-8,mt+i*ch+15,{"font-size":11,"text-anchor":"end"});cols.forEach((col,j)=>{const cell=value(row,col)||{};const node=S("rect",{x:ml+j*cw,y:mt+i*ch,width:cw-2,height:ch-2,fill:cell.color||"#edf1f5",stroke:cell.state==="NOT_APPLICABLE"?"#8f9caa":"white","stroke-dasharray":cell.state==="NOT_APPLICABLE"?"2 2":"none"});interactive(node,cell.tip||`${esc(row)} · ${esc(col)}: ${fmt(cell.v)}`,opts.onClick?()=>opts.onClick(row,col):null);svg.append(node);});});svg._reset=zoomable(svg,w,h);if(rows.length<rowsL.length||cols.length<colsL.length)svg._note=`Showing ${rows.length}/${rowsL.length} rows and ${cols.length}/${colsL.length} columns; filter or use table/exports for all data.`;return svg;
}
function boxPlot(groups,label="Observed value distribution") {
  const items=[];
  for(const [name,values] of groups){
    const sorted=values.map(number).filter(v=>v!==null).sort((a,b)=>a-b);if(!sorted.length)continue;
    const q=p=>{const at=(sorted.length-1)*p,lo=Math.floor(at),hi=Math.ceil(at);return sorted[lo]+(sorted[hi]-sorted[lo])*(at-lo);};
    items.push({name,min:sorted[0],q1:q(.25),median:q(.5),q3:q(.75),max:sorted[sorted.length-1],n:sorted.length});
  }
  if(!items.length)return empty("No supplied point estimates are available for this distribution.");
  const shown=items.slice(0,40),w=850,h=75+shown.length*42,left=250,right=100;
  const lo=Math.min(...items.map(r=>r.min)),hi=Math.max(...items.map(r=>r.max)),x=v=>left+(v-lo)/(hi-lo||1)*(w-left-right);
  const svg=S("svg",{viewBox:`0 0 ${w} ${h}`,width:w,role:"img","aria-label":label+" box plot"});
  for(let k=0;k<=4;k++){const value=lo+(hi-lo)*k/4;axisText(svg,fmt(value,2),x(value),h-28,{"text-anchor":"middle"});}
  shown.forEach((r,i)=>{const y=25+i*42;axisText(svg,r.name,left-12,y+4,{"text-anchor":"end"});
    svg.append(S("line",{x1:x(r.min),x2:x(r.max),y1:y,y2:y,stroke:"#276685"}));
    for(const value of [r.min,r.max])svg.append(S("line",{x1:x(value),x2:x(value),y1:y-7,y2:y+7,stroke:"#276685"}));
    const box=S("rect",{x:x(r.q1),y:y-10,width:Math.max(2,x(r.q3)-x(r.q1)),height:20,fill:"#bcd7e7",stroke:"#276685"});
    interactive(box,`${esc(r.name)}; n=${r.n}; median ${fmt(r.median)}; Q1–Q3 ${fmt(r.q1)}–${fmt(r.q3)}; min–max ${fmt(r.min)}–${fmt(r.max)}. Whiskers show observed min/max.`);
    svg.append(box,S("line",{x1:x(r.median),x2:x(r.median),y1:y-10,y2:y+10,stroke:"#14324a","stroke-width":2}));axisText(svg,"n="+r.n,w-right+15,y+4);
  });
  axisText(svg,label,w/2,h-5,{"text-anchor":"middle"});svg._reset=zoomable(svg,w,h);
  svg._note="Supplied point estimates: box = Q1–Q3; line = median; whiskers = observed min/max. These are descriptive, not confidence or posterior intervals."+(items.length>40?" First 40 groups drawn; exports retain all records.":"");return svg;
}
function lineChart(rows,{xField="month",yField="prevalence",groupField="country",label="Monthly prevalence"}={}) {
  const months=unique(rows.map(r=>r[xField])).sort();if(!months.length)return empty("No recorded time periods are available.");const w=800,h=400,m={l:65,r:140,t:25,b:85};const valid=rows.filter(r=>number(r[yField])!==null);if(!valid.length)return empty("No supplied prevalence estimates meet the upstream reporting criteria; missing months are not zero.");const max=Math.max(1,...valid.map(r=>number(r[yField]))),sx=month=>m.l+months.indexOf(month)/Math.max(1,months.length-1)*(w-m.l-m.r),sy=value=>h-m.b-value/max*(h-m.t-m.b);const svg=S("svg",{viewBox:`0 0 ${w} ${h}`,width:w,role:"img","aria-label":label});const colors=["#24759e","#9e5a39","#4f8a67","#785e9a"];for(let k=0;k<=4;k++){const val=max*k/4;svg.append(S("line",{x1:m.l,x2:w-m.r,y1:sy(val),y2:sy(val),stroke:"#e2eaf0"}));axisText(svg,fmt(val,2),m.l-8,sy(val)+4,{"text-anchor":"end"});}months.forEach(month=>axisText(svg,month,sx(month),h-m.b+18,{"font-size":10,transform:`rotate(-35 ${sx(month)} ${h-m.b+18})`,"text-anchor":"end"}));let gi=0;for(const [group,list]of groupBy(rows,groupField)){const color=colors[gi%colors.length];let previous=null;for(const month of months){const record=list.find(r=>r[xField]===month);const value=record?number(record[yField]):null;if(value===null){previous=null;continue;}const point={x:sx(month),y:sy(value)};if(previous)svg.append(S("line",{x1:previous.x,y1:previous.y,x2:point.x,y2:point.y,stroke:color,"stroke-width":2}));const circle=S("circle",{cx:point.x,cy:point.y,r:5,fill:color});interactive(circle,`${esc(group)} · ${esc(month)}; prevalence ${fmt(value)}; numerator ${fmt(record.mutant_or_mixed)}; callable denominator ${fmt(record.callable)}; samples ${fmt(record.samples)}`);svg.append(circle);previous=point;}axisText(svg,group,w-m.r+12,30+gi*18,{fill:color});gi++;}axisText(svg,"Recorded month",w/2,h-4,{"text-anchor":"middle"});axisText(svg,"Prevalence",13,h/2,{transform:`rotate(-90 13 ${h/2})`});svg._note="Gaps represent missing/insufficient upstream estimates; no interpolation or zero substitution.";return svg;
}
function numericColor(value,max=1) { const n=number(value);return n===null?"#edf1f5":`rgba(31,105,148,${Math.max(.12,Math.min(1,.12+.88*n/(max||1)))})`; }
const STATUS_COLOR = { PASS: "#0ca30c", WARNING: "#fab219", FAIL: "#d03b3b", NOT_EVALUATED: "#c3cad3", NOT_APPLICABLE: "#eceff3", UNKNOWN: "#8795ac",
  MUTANT: "#d03b3b", MIXED: "#fab219", WILD_TYPE: "#0ca30c", BELOW_REPORTING_FLOOR: "#9fb4c7",
  AMPLIFIED: "#d03b3b", DELETED: "#d03b3b", SINGLE_COPY: "#0ca30c", PRESENT: "#0ca30c",
  // Clonality is a category, not a quality state: categorical slots 1-2 (blue / orange, validated
  // for colour-vision deficiency), gray when undecided. Charts always label them.
  MONOCLONAL: "#2a78d6", POLYCLONAL: "#eb6834", INDETERMINATE: "#a7b1bd" };

/* ---------- map (equirectangular, bundled Natural Earth) ---------- */
function mapView(points, opts = {}) { points = points.filter(p => number(p.lat) !== null && number(p.lon) !== null && Math.abs(number(p.lat)) <= 90 && Math.abs(number(p.lon)) <= 180); if (!points.length) return empty("Map unavailable: no valid recorded geographic coordinates for the selected data."); if (!D.geo || !Array.isArray(D.geo.features)) return empty("Map unavailable: bundled geographic boundaries were not supplied; no schematic boundaries substituted."); const w = 760, h = 420; const svg = S("svg", { viewBox: `0 0 ${w} ${h}`, width: w, role: "img", "aria-label": "map" });
  const lons = points.map(p => p.lon), lats = points.map(p => p.lat); const pad = 8;
  let [lo0, lo1, la0, la1] = points.length ? [Math.min(...lons) - pad, Math.max(...lons) + pad, Math.min(...lats) - pad, Math.max(...lats) + pad] : [-20, 55, -36, 38];
  const px = lon => (lon - lo0) / (lo1 - lo0) * w, py = lat => (la1 - lat) / (la1 - la0) * h;
  if (D.geo) D.geo.features.forEach(f => { const g = f.geometry; const polys = g.type === "MultiPolygon" ? g.coordinates : [g.coordinates];
    polys.forEach(poly => { const d = poly.map(ring => "M" + ring.map(([x, y]) => px(x).toFixed(1) + "," + py(y).toFixed(1)).join("L") + "Z").join(""); const p = S("path", { d, fill: "#eef2f6", stroke: "#9aa8b8", "stroke-width": 0.6 }); p.addEventListener("mousemove", e => showTip(e, esc(f.properties.name))); p.addEventListener("mouseleave", hideTip); svg.append(p); }); });
  points.forEach(p => { const c = S("circle", { cx: px(p.lon), cy: py(p.lat), r: p.r || 6, fill: p.color || "#c62828", "fill-opacity": 0.8, stroke: "#fff", tabindex: 0 }); c.addEventListener("mousemove", e => showTip(e, p.tip)); c.addEventListener("mouseleave", hideTip); if (opts.onClick) c.addEventListener("click", () => opts.onClick(p)); svg.append(c); });
  const note = S("text", { x: 6, y: h - 6, "font-size": 10, fill: "#5b6575" }); note.textContent = "Natural Earth 1:110m (bundled, offline); positions aggregated"; svg.append(note);
  svg._reset = zoomable(svg, w, h); return svg; }

/* ---------- tree (Newick, rectangular, collapsible) ---------- */
function parseNewick(s) { let i = 0; function node() { const n = { children: [], name: "", len: null }; if (s[i] === "(") { i++; do { n.children.push(node()); } while (s[i] === "," && ++i); i++; }
    let lab = ""; while (i < s.length && !",():;".includes(s[i])) lab += s[i++]; n.name = lab; if (s[i] === ":") { i++; let l = ""; while (i < s.length && !",();".includes(s[i])) l += s[i++]; n.len = number(l); } return n; } return node(); }
function treeView(newick, meta = {}) { const root = parseNewick(newick.trim()); const w = 760; let leaves = 0; const svg = S("svg", { role: "img", width: w, "aria-label": "phylogenetic tree" });
  const allNodes = n => [n, ...n.children.flatMap(allNodes)]; const topologyOnly = allNodes(root).some(n => n !== root && n.len === null); svg._note = topologyOnly ? "Topological layout only: some branch lengths were not supplied; horizontal spacing is not evolutionary distance." : "Horizontal branch lengths follow the supplied Newick distances; internal labels retain supplied support values."; function draw() { svg.replaceChildren(); leaves = 0; const depth = (n, d) => { n.d = d + (topologyOnly ? 1 : (n.len ?? 0)); n.children.forEach(c => depth(c, n.d)); }; depth(root, 0);
    const maxd = Math.max(...(function all(n) { return [n.d, ...n.children.flatMap(all)]; })(root)) || 1;
    const place = (n) => { if (!n.children.length || n.collapsed) { n.y = 20 + leaves++ * 22; } else { n.children.forEach(place); n.y = (n.children[0].y + n.children[n.children.length - 1].y) / 2; } }; place(root);
    const h = 40 + leaves * 22; svg.setAttribute("viewBox", `0 0 ${w} ${h}`); const sx = d => 20 + d / maxd * (w - 220);
    const rec = (n) => { const x = sx(n.d); if (n.children.length && !n.collapsed) { svg.append(S("line", { x1: x, x2: x, y1: n.children[0].y, y2: n.children[n.children.length - 1].y, stroke: "#333" }));
        n.children.forEach(c => { svg.append(S("line", { x1: x, x2: sx(c.d), y1: c.y, y2: c.y, stroke: "#333" })); rec(c); });
        if (n.name) { const t = S("text", { x: x + 3, y: n.y - 3, "font-size": 9, fill: "#5b6575" }); t.textContent = n.name; svg.append(t); }
        const k = S("circle", { cx: x, cy: n.y, r: 4, fill: "#1f5f8b", tabindex: 0 }); k.addEventListener("click", () => { n.collapsed = true; draw(); }); k.addEventListener("mousemove", e => showTip(e, "support: " + esc(n.name || "—") + "<br>click to collapse")); k.addEventListener("mouseleave", hideTip); svg.append(k); }
      else { const lab = n.collapsed ? `[collapsed clade: ${(function c(x) { return x.children.length ? x.children.reduce((a, b) => a + c(b), 0) : 1; })(n)} tips]` : n.name;
        const t = S("text", { x: x + 4, y: n.y + 4, "font-size": 12, "font-weight": STATE.selectedSample === n.name ? "700" : "400", tabindex: 0 }); t.textContent = lab + (meta[n.name] ? "  · " + meta[n.name] : "");
        t.addEventListener("click", () => { if (n.collapsed) { n.collapsed = false; draw(); } else selectSample(n.name); }); svg.append(t); } };
    rec(root); svg._reset = zoomable(svg, w, h); }
  draw(); return svg; }

/* ---------- network (circular layout) ---------- */
function networkView(nodes, edges, threshold) { const original = nodes.length; nodes = nodes.slice(0,60); edges = edges.filter(e => nodes.includes(e.a) && nodes.includes(e.b) && number(e.w) !== null).slice(0,600); const w = 560, h = 460, cx = w / 2, cy = h / 2, R = 170; const svg = S("svg", { viewBox: `0 0 ${w} ${h}`, width: w, role: "img", "aria-label": "relatedness network" });
  const pos = {}; nodes.forEach((n, i) => pos[n] = [cx + R * Math.cos(2 * Math.PI * i / nodes.length), cy + R * Math.sin(2 * Math.PI * i / nodes.length)]);
  edges.filter(e => e.w >= threshold).forEach(e => { const l = S("line", { x1: pos[e.a][0], y1: pos[e.a][1], x2: pos[e.b][0], y2: pos[e.b][1], stroke: "#1f5f8b", "stroke-width": Math.max(.5, 1 + 6 * e.w), "stroke-opacity": 0.6 }); l.addEventListener("mousemove", ev => showTip(ev, `${esc(e.a)}–${esc(e.b)}: ${fmt(e.w)} (${esc(e.method)})`)); l.addEventListener("mouseleave", hideTip); svg.append(l); });
  nodes.forEach(n => { const c = S("circle", { cx: pos[n][0], cy: pos[n][1], r: 7, fill: STATE.selectedSample === n ? "#c62828" : "#14324a", tabindex: 0 }); c.addEventListener("click", () => selectSample(n)); svg.append(c); const t = S("text", { x: pos[n][0] + 9, y: pos[n][1] + 4, "font-size": 11 }); t.textContent = n; svg.append(t); });
  svg._reset = zoomable(svg, w, h); if (nodes.length < original) svg._note = `First ${nodes.length} of ${original} nodes shown; filters and table retain the full data.`; return svg; }

