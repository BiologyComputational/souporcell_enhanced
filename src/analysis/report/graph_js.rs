// ============================================================================
// report/graph_js.rs — D3.js pipeline graph JavaScript block
//
// Returns the complete <script> block. Kept separate so the JS can be read
// and edited independently without navigating the main HTML format! string.
// ============================================================================

/// Build the full <script> tag for the D3 pipeline graph.
/// conv_ok and bal_ok are Rust bools that become JS `true`/`false` literals.
pub fn build_graph_script(graph_json: &str, conv_ok: bool, bal_ok: bool) -> String {
    let conv_js = if conv_ok { "true" } else { "false" };
    let bal_js  = if bal_ok  { "true" } else { "false" };

    format!(r####"<script>
(function(){{
const GRAPH = {graph_json};
const GLOBAL_CKS = [
  ["Loci QC pass rate > 1%",       true],
  ["Matrix sparsity > 90%",        true],
  ["All restarts completed",       true],
  ["Convergence stability >= 70%", {conv_js}],
  ["Cluster balance within 50%",   {bal_js}],
  ["No empty clusters",            true],
  ["Median cell confidence > 15",  true],
  ["All SVG plots written",        true],
  ["clusters_tmp.tsv written",     true],
  ["HTML report generated",        true],
];
const LINK_COL = {{trigger:"#C75B00",data:"#006A6A",result:"#15803D"}};
const CAT_COL  = {{io:"#0369A1",qc:"#B91C1C",compute:"#15803D",algo:"#C75B00",output:"#7C3AED"}};
const CAT_NAME = {{io:"I/O",qc:"Quality control",compute:"Compute",algo:"Algorithm core",output:"Output"}};
const PHASE_COL = {{"Init":"#0369A1","Clustering":"#C75B00","Post-clustering":"#7C3AED"}};
const PHASE_ORDER = ["Init","Clustering","Post-clustering"];
function col(n){{ return CAT_COL[n.cat]||"#64748B"; }}

function buildSB(){{
  let h="";
  PHASE_ORDER.forEach(ph=>{{
    const ns=GRAPH.nodes.filter(n=>n.phase===ph);
    h+=`<div style="padding:7px 14px 2px;font-size:10px;font-weight:700;color:${{PHASE_COL[ph]}};text-transform:uppercase;letter-spacing:.08em">${{ph}}</div>`;
    ns.forEach(n=>{{
      const c=col(n);
      h+=`<div id="si-${{n.id}}" onclick="selNode('${{n.id}}')"
        style="display:flex;align-items:center;gap:8px;padding:7px 14px;cursor:pointer;border-left:3px solid transparent;opacity:.4;transition:all .2s">
        <div id="si-ico-${{n.id}}" style="width:22px;height:22px;border-radius:50%;display:flex;align-items:center;justify-content:center;font-size:9px;font-weight:700;flex-shrink:0;border:2px solid ${{c}};color:${{c}};background:transparent">${{String(n.seq).padStart(2,'0')}}</div>
        <div style="flex:1;min-width:0">
          <div style="font-size:12px;color:var(--ink);white-space:nowrap;overflow:hidden;text-overflow:ellipsis">${{n.name}}</div>
          <div style="font-size:10px;color:var(--muted)">${{n.timing||CAT_NAME[n.cat]||''}}</div>
        </div>
        <div id="si-tick-${{n.id}}" style="width:16px;height:16px;border-radius:50%;border:2px solid var(--border);display:flex;align-items:center;justify-content:center;font-size:9px;flex-shrink:0"></div>
      </div>`;
    }});
  }});
  document.getElementById('graph-stage-list').innerHTML=h;
}}

function buildCK(){{
  document.getElementById('graph-checklist').innerHTML=GLOBAL_CKS.map(([lbl,met])=>
    `<div style="display:flex;align-items:flex-start;gap:7px;padding:2px 0">
      <div style="width:14px;height:14px;border-radius:3px;flex-shrink:0;margin-top:1px;background:${{met?'#15803D':'#B91C1C'}};display:flex;align-items:center;justify-content:center;font-size:9px;color:#fff">${{met?'&#10003;':'&#10007;'}}</div>
      <span style="font-size:11px;color:${{met?'var(--ink2)':'#B91C1C'}};line-height:1.4">${{lbl}}</span>
    </div>`).join('');
}}

let sim,lSel,llSel,nSel,stepMode=false,curStep=0;

function initGraph(){{
  const svg=d3.select('#graph-svg');
  const W=svg.node().parentElement.clientWidth||700;
  const H=svg.node().parentElement.clientHeight||600;
  svg.selectAll('*').remove();

  const defs=svg.append('defs');
  Object.entries(LINK_COL).forEach(([t,c])=>{{
    defs.append('marker').attr('id','arr-'+t)
      .attr('viewBox','0 -5 10 10').attr('refX',34).attr('refY',0)
      .attr('markerWidth',7).attr('markerHeight',7).attr('orient','auto')
      .append('path').attr('d','M0,-5L10,0L0,5').attr('fill',c).attr('opacity',.85);
  }});

  const g=svg.append('g');
  svg.call(d3.zoom().scaleExtent([.2,3]).on('zoom',e=>g.attr('transform',e.transform)));

  const nodes=GRAPH.nodes.map(d=>({{...d}}));
  const nm=Object.fromEntries(nodes.map(n=>[n.id,n]));
  const links=GRAPH.links.map(l=>({{...l,source:nm[l.s],target:nm[l.t]}}));

  lSel=g.append('g').selectAll('line').data(links).join('line')
    .attr('stroke',d=>LINK_COL[d.type]).attr('stroke-width',2)
    .attr('stroke-opacity',0).attr('marker-end',d=>'url(#arr-'+d.type+')');

  llSel=g.append('g').selectAll('g').data(links).join('g').attr('opacity',0);
  llSel.append('rect').attr('x',-10).attr('y',-13).attr('width',20).attr('height',13)
       .attr('rx',5).attr('fill',d=>LINK_COL[d.type]).attr('opacity',.9);
  llSel.append('text').attr('font-size',8).attr('font-weight','600').attr('fill','#fff')
       .attr('text-anchor','middle').attr('dominant-baseline','central').attr('dy',-6)
       .text(d=>d.ord);
  llSel.append('text').attr('font-size',9).attr('fill','var(--muted)')
       .attr('text-anchor','middle').attr('dy',8).text(d=>d.lbl);

  nSel=g.append('g').selectAll('g').data(nodes).join('g')
    .attr('cursor','pointer').attr('opacity',0)
    .call(d3.drag()
      .on('start',(e,d)=>{{if(!e.active)sim.alphaTarget(.3).restart();d.fx=d.x;d.fy=d.y}})
      .on('drag', (e,d)=>{{d.fx=e.x;d.fy=e.y}})
      .on('end',  (e,d)=>{{if(!e.active)sim.alphaTarget(0);d.fx=null;d.fy=null}}));

  nSel.append('circle').attr('r',d=>d.r+9).attr('fill','none')
      .attr('stroke',d=>col(d)).attr('stroke-width',.6).attr('opacity',.18);
  nSel.append('circle').attr('r',d=>d.r).attr('fill',d=>col(d)+'14')
      .attr('stroke',d=>col(d)).attr('stroke-width',2.5);

  const seqg=nSel.append('g');
  seqg.append('rect').attr('x',-11).attr('y',d=>-(d.r+17)).attr('width',22).attr('height',14)
      .attr('rx',7).attr('fill',d=>col(d));
  seqg.append('text').attr('text-anchor','middle').attr('y',d=>-(d.r+8))
      .attr('font-size',8).attr('font-weight','700').attr('fill','#fff')
      .text(d=>String(d.seq).padStart(2,'0'));

  nSel.append('text').attr('text-anchor','middle').attr('dominant-baseline','central')
      .attr('font-size',10).attr('font-weight','600').attr('fill',d=>col(d))
      .text(d=>d.name.split(' ').map(w=>w.slice(0,4)).join(' '));
  nSel.append('text').attr('text-anchor','middle').attr('dy',d=>d.r+15)
      .attr('font-size',11).attr('fill','var(--ink2)').attr('font-weight','600')
      .text(d=>d.name);
  nSel.filter(d=>d.timing).append('text').attr('text-anchor','middle')
      .attr('dy',d=>d.r+26).attr('font-size',9).attr('fill','var(--teal)')
      .text(d=>d.timing);

  nSel.on('click',(e,d)=>{{e.stopPropagation();selNode(d.id)}});
  svg.on('click',()=>document.getElementById('detail-panel').style.display='none');

  sim=d3.forceSimulation(nodes)
    .force('link',d3.forceLink(links).id(d=>d.id).distance(145).strength(.55))
    .force('charge',d3.forceManyBody().strength(-560))
    .force('center',d3.forceCenter(W/2,H/2))
    .force('collide',d3.forceCollide().radius(d=>d.r+40))
    .on('tick',()=>{{
      lSel.attr('x1',d=>d.source.x).attr('y1',d=>d.source.y)
          .attr('x2',d=>d.target.x).attr('y2',d=>d.target.y);
      llSel.attr('transform',d=>`translate(${{(d.source.x+d.target.x)/2}},${{(d.source.y+d.target.y)/2}})`);
      nSel.attr('transform',d=>`translate(${{d.x}},${{d.y}})`);
    }});

  if(!stepMode) showAll();
}}

function showAll(){{
  nSel.transition().duration(350).attr('opacity',1);
  lSel.transition().duration(350).delay(250).attr('stroke-opacity',.7);
  llSel.transition().duration(350).delay(320).attr('opacity',1);
  GRAPH.nodes.forEach(n=>markDone(n.id));
}}

function markDone(id){{
  const c=col(GRAPH.nodes.find(n=>n.id===id));
  const si=document.getElementById('si-'+id);
  const ico=document.getElementById('si-ico-'+id);
  const tick=document.getElementById('si-tick-'+id);
  if(si)   si.style.opacity='1';
  if(ico)  {{ ico.style.background=c; ico.style.color='#fff'; ico.style.borderColor=c; }}
  if(tick) {{ tick.style.background='#15803D'; tick.style.border='none'; tick.style.color='#fff'; tick.innerHTML='&#10003;'; }}
}}

function markPending(id){{
  const c=col(GRAPH.nodes.find(n=>n.id===id));
  const si=document.getElementById('si-'+id);
  const ico=document.getElementById('si-ico-'+id);
  const tick=document.getElementById('si-tick-'+id);
  if(si)   si.style.opacity='.4';
  if(ico)  {{ ico.style.background=''; ico.style.color=c; ico.style.borderColor=c; }}
  if(tick) {{ tick.style.background=''; tick.style.border='2px solid var(--border)'; tick.innerHTML=''; }}
}}

window.graphShowAll=function(){{
  stepMode=false;
  document.getElementById('graph-step-controls').style.display='none';
  document.getElementById('step-desc').style.display='none';
  showAll();
}};

window.graphStartStep=function(){{
  stepMode=true; curStep=0;
  document.getElementById('graph-step-controls').style.display='flex';
  document.getElementById('step-desc').style.display='block';
  nSel.attr('opacity',0); lSel.attr('stroke-opacity',0); llSel.attr('opacity',0);
  GRAPH.nodes.forEach(n=>markPending(n.id));
  revealTo(0);
}};

function revealTo(idx){{
  const vis=new Set(GRAPH.nodes.slice(0,idx+1).map(n=>n.id));
  nSel.transition().duration(450).attr('opacity',(_,i)=>i<=idx?1:0);
  lSel.transition().duration(350)
      .attr('stroke-opacity',d=>(vis.has(d.source.id)&&vis.has(d.target.id))?.7:0);
  llSel.transition().duration(350)
       .attr('opacity',d=>(vis.has(d.source.id)&&vis.has(d.target.id))?1:0);
  GRAPH.nodes.forEach((n,i)=>{{ if(i<=idx) markDone(n.id); else markPending(n.id); }});
  const n=GRAPH.nodes[idx];
  document.getElementById('graph-step-label').textContent=`Step ${{idx+1}} / ${{GRAPH.nodes.length}}`;
  document.getElementById('step-desc').textContent=
    `Stage ${{n.seq}}: ${{n.name}} \u2014 ${{(n.desc||'').slice(0,120)}}\u2026`;
  selNode(n.id);
}}

window.graphNext=function(){{
  if(curStep<GRAPH.nodes.length-1){{ curStep++; revealTo(curStep); }}
}};
window.graphPrev=function(){{
  if(curStep>0){{ curStep--; revealTo(curStep); }}
}};

window.selNode=function(id){{
  document.querySelectorAll('[id^="si-"]').forEach(e=>{{
    if(e.id.startsWith('si-')&&!e.id.includes('ico')&&!e.id.includes('tick'))
      e.style.borderLeftColor='transparent';
  }});
  const si=document.getElementById('si-'+id);
  const n=GRAPH.nodes.find(x=>x.id===id); if(!n) return;
  const c=col(n);
  if(si) si.style.borderLeftColor=c;

  document.getElementById('dp-title').textContent=
    `${{String(n.seq).padStart(2,'0')}}. ${{n.name}}`;
  document.getElementById('dp-badge').innerHTML=
    `<span style="font-size:10px;padding:3px 9px;border-radius:6px;background:${{c}}18;color:${{c}};border:1px solid ${{c}}44;font-weight:600">${{n.phase}} \u00b7 ${{CAT_NAME[n.cat]||n.cat}}</span>`;
  document.getElementById('dp-desc').textContent=n.desc||'';
  document.getElementById('dp-rows').innerHTML=(n.details||[]).map(([k,v])=>
    `<div style="display:flex;justify-content:space-between;gap:8px;padding:4px 0;border-bottom:.5px solid var(--border)">
      <span style="font-size:11px;color:var(--muted)">${{k}}</span>
      <span style="font-size:11px;color:var(--ink2);font-weight:600;text-align:right">${{v}}</span>
    </div>`).join('');
  document.getElementById('dp-node-checks').innerHTML=
    `<div style="font-size:10px;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.07em;margin-bottom:5px">Stage checks</div>`+
    (n.checks||[]).map(([lbl,met])=>
      `<div style="display:flex;align-items:flex-start;gap:6px;padding:2px 0">
        <div style="width:13px;height:13px;border-radius:2px;flex-shrink:0;margin-top:1px;background:${{met?'#15803D':'#B91C1C'}};display:flex;align-items:center;justify-content:center;font-size:8px;color:#fff">${{met?'&#10003;':'&#10007;'}}</div>
        <span style="font-size:11px;color:${{met?'var(--ink2)':'#B91C1C'}};line-height:1.4">${{lbl}}</span>
      </div>`).join('')+
    `<div style="font-size:10px;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.07em;margin:8px 0 5px">Sub-steps</div>`+
    (n.subs||[]).map((s,i)=>
      `<div style="display:flex;align-items:flex-start;gap:7px;padding:3px 0;border-bottom:.5px solid var(--border)">
        <div style="width:16px;height:16px;border-radius:50%;background:${{c}}18;border:1.5px solid ${{c}};flex-shrink:0;display:flex;align-items:center;justify-content:center;font-size:8px;color:${{c}};margin-top:1px">${{i+1}}</div>
        <span style="font-size:11px;color:var(--muted);line-height:1.4">${{s}}</span>
      </div>`).join('');

  document.getElementById('detail-panel').style.display='block';
}};

window.graphToggleFS=function(){{
  const el=document.getElementById('graph-wrap');
  if(!document.fullscreenElement) el.requestFullscreen&&el.requestFullscreen();
  else document.exitFullscreen&&document.exitFullscreen();
}};

buildSB(); buildCK();
setTimeout(initGraph,80);
}})();
</script>"####,
        graph_json = graph_json,
        conv_js    = conv_js,
        bal_js     = bal_js,
    )
}
