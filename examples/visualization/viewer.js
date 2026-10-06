"use strict";

// The scientific judgments and coordinates come from the embedded saved-data export.
const PAYLOAD = JSON.parse(document.getElementById("example-payload").textContent);
const EXAMPLES = PAYLOAD.data.examples;
const COLORS = {H:"#f4f5f7",C:"#424b58",N:"#3564cb",O:"#d04a4a",F:"#51a568",B:"#d89575"};
const VIEWERS = [];
let selectedExample = null;
let profileViewer = null;
let currentProfile = null;

const $ = id => document.getElementById(id);
const esc = value => String(value ?? "").replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const number = (value, digits=2) => Number.isFinite(Number(value)) && value !== null && value !== undefined ? Number(value).toFixed(digits) : "未算出";
const text = value => typeof value === "string" ? value : JSON.stringify(value ?? "未記録");
const table = (headers, rows) => `<div class="table-wrap"><table><thead><tr>${headers.map(h=>`<th>${esc(h)}</th>`).join("")}</tr></thead><tbody>${rows.map(row=>`<tr>${row.map(v=>`<td>${esc(v)}</td>`).join("")}</tr>`).join("")}</tbody></table></div>`;

function caseKind(ex) {
  const id = (ex.id + " " + ex.system).toLowerCase();
  if (id.includes("tma")) return "tma";
  if (id.includes("bh3")) return "bh3";
  if (id.includes("oh") && id.includes("ch4")) return "oh";
  if (id.includes("hf")) return "exchange";
  return "inversion";
}

const COPY = {
  inversion: {
    short:"NH₃反転", title:"NH₃の反転 — 同じ状態を結ぶTSを認証する",
    purpose:"平面の入力は極小ではありません。虚振動に沿った両側の降下と停留点の検証によって、ピラミッドの反転を確認します。",
    outcome:"縮退反転を認証",
    flow:[["平面の構造","入力の停留性を確認"],["虚振動1本","反応方向に沿って下る"],["両側の接続","上下の構造を対応付ける"],["TSを認証","同じ化学状態の反転"]],
    conclusion:"上下のピラミッドは同じ化学状態に属します。その間のTSと両側の接続があるため、同じ状態に戻るという理由だけで無反応とは判断しません。",
    limits:"縮退した反転であり、別の生成物を発見した例ではありません。保存された反対側のmode-follow構造を3Dで表示します。",
    structureNote:"保存された上下の極小とTS。原子の元素・個数・番号を保持しています。"
  },
  exchange: {
    short:"NH₃・HF", title:"NH₃・HF — 探索から水素交換を見つける",
    purpose:"生成物を指定しないpilotから、NH₃・HF内で水素を交換するTSが得られました。探索候補と認証した反応を分けて追います。",
    outcome:"水素交換を認証",
    flow:[["候補を探索","結合編集とNT2"],["TSを精密化","DFTの停留点へ"],["振動・QRC","両側と原子対応を確認"],["条件付き比較","錯体を基準にGを評価"]],
    conclusion:"化学式と化学状態は同じでも、原子対応を追うと水素交換が分かります。TS、虚振動、端点接続をまとめて認証しています。",
    limits:"TMA・HFを含むpilotのうち、この認証経路はNH₃・HF由来です。TMAの反応を認証した結果として扱いません。W4のQM結果をW5の熱化学で再評価しています。",
    structureNote:"QRCの一方は独立した最適化、反対側は実装が再利用した対称像です。対称像には出所を表示し、独立したQM計算とは区別します。"
  },
  tma: {
    short:"TMA・(HF)₂", title:"TMA・(HF)₂ — 異なる入力が同じ極小へ合流する",
    purpose:"neutralとshared_protonという別名で与えた構造を、保存されたDFT極小と対比します。ラベルと独立した極小の違いを示す例です。",
    outcome:"同じbasinに合流",
    flow:[["2種類の入力","宣言した端点候補"],["DFT最適化","極小を検証する"],["同一性を照合","同じbasinへ合流"],["順位対象なし","認証した別TSがない"]],
    conclusion:"別名の入力だけでは別の反応物・生成物とは言えません。この保存結果の宣言反応はsame_basinであり、TSと障壁の順位を持ちません。",
    limits:"探索予算内で認証した別経路がないことは、あらゆる反応の不存在の証明ではありません。表示は完了したW4結果。W5 freshは未完了です。",
    structureNote:"入力候補と到達した極小を区別して表示します。入力構造に自由エネルギーの極小値を与えたり、架空のTSを置いたりしません。"
  },
  bh3: {
    short:"BH₃＋NH₃", title:"BH₃＋NH₃ — 会合を障壁なしの経路として扱う",
    purpose:"ルイス酸・塩基の会合について、保存された走査節点と中点追加を確認します。会合の自由エネルギーと経路上の電子エネルギーを分けます。",
    outcome:"解像度内で障壁なし",
    flow:[["分離の単量体","共通のゼロを置く"],["会合経路を走査","電子エネルギーの節点"],["中点を追加","節点間の山を確認"],["捕獲律速として記録","通常の障壁順位から分ける"]],
    conclusion:"調べた計算面・節点の分解能内で山が検出されず、barrierless_at_resolutionと判断しています。通常の一次鞍点を経由する反応とは表示を分けます。",
    limits:"速度定数は未計算です。反応の自由エネルギーが負であることと、実際の速度や生成物比が決まることは別です。",
    structureNote:"最適化したBH₃、NH₃、会合体を表示します。宣言した接近配置も選択できます。自由エネルギー図のRは分離した単量体の和で、接近配置の構造とは区別します。"
  },
  oh: {
    short:"OH＋CH₄", title:"OH＋CH₄ — 複数の候補を個別反応として検証する",
    purpose:"ラジカルによる水素引き抜きを中心に、探索、TSと接続、電子状態、分離基準のエネルギーを追います。",
    outcome:"引き抜き経路を認証",
    flow:[["組成から探索","配置とNT2候補"],["DFT入口を選ぶ","低レベルの候補を検証"],["TS・QRC・電子状態","採用と未解決を分ける"],["個別反応を比較","出発状態の違いを残す"]],
    conclusion:"引き抜きのTSからCH₃＋H₂O側へ接続します。錯体のGと分離したOH＋CH₄のGを区別し、同じゼロから反応量を読むことができます。",
    limits:"ほかの保存反応は出発状態が異なります。順位をそのままOH＋CH₄からの主要生成物の順位へ読み替えません。W5で更新したのは熱化学です。",
    structureNote:"主表示は認証した水素引き抜きの反応前・TS・反応後。原子対応を保ち、水素の位置とO–H/C–H距離を比較できます。"
  }
};

function xyzOf(structure) {
  if (structure.xyz) return structure.xyz;
  return `${structure.atoms.length}\n${structure.label || structure.role}\n` + structure.atoms.map(a=>`${a.element} ${a.x} ${a.y} ${a.z}`).join("\n") + "\n";
}

function validStructures(ex) {
  return (ex.structures || []).filter(s=>s.atoms && s.atoms.length);
}

function makeViewer(element) {
  if (!window.$3Dmol) throw new Error("3Dmol.jsが読み込めませんでした。");
  return $3Dmol.createViewer(element,{backgroundColor:"white",antialias:true});
}

function putStructure(viewer, structure, distanceId) {
  viewer.removeAllModels(); viewer.removeAllLabels(); viewer.removeAllShapes();
  $(distanceId).textContent="原子を2つクリックすると距離を表示します。";
  const model = viewer.addModel(xyzOf(structure), "xyz");
  viewer.setStyle({}, {stick:{radius:0.12},sphere:{scale:0.27}});
  for (const [elem,color] of Object.entries(COLORS)) viewer.setStyle({elem},{stick:{radius:0.12,color},sphere:{scale:0.27,color}});
  if ($("atom-labels").checked) {
    for (const atom of model.selectedAtoms({})) viewer.addLabel(`${atom.elem}${atom.index+1}`,{position:atom,fontSize:12,fontColor:"#233e46",backgroundOpacity:0,inFront:true});
  }
  let previous = null;
  viewer.setClickable({},true,atom=>{
    if (previous !== null && previous.index !== atom.index) {
      const d=Math.hypot(atom.x-previous.x,atom.y-previous.y,atom.z-previous.z);
      $(distanceId).textContent=`${previous.elem}${previous.index+1} – ${atom.elem}${atom.index+1}: ${d.toFixed(3)} Å`;
      viewer.removeAllShapes();
      viewer.addLine({start:previous,end:atom,color:"#176f75",dashed:true,linewidth:2});
      previous=null;
    } else {
      previous=atom;
      $(distanceId).textContent=`${atom.elem}${atom.index+1} を選択。比較する原子をクリックしてください。`;
    }
    viewer.render();
  });
  viewer.zoomTo();
  viewer.zoom(structure.atoms.length<10?1.8:1.35);
  viewer.rotate(65,{x:1,y:0,z:0});
  viewer.render();
  return viewer;
}

function structureInfo(s) {
  let note = s.observed_role || s.role || "保存された構造";
  const normalized=note.toLowerCase();
  if (/symmetry|image/.test(normalized)) note="計算結果の対称像（独立QMではありません）";
  else if (/input|seed|declared/.test(normalized)) note="入力候補（認証した極小・TSとは区別）";
  else if (/saddle|\bts\b/.test(normalized)) note="検証したTS";
  else if (/minimum|basin|optimized|qrc/.test(normalized)) note="保存された最適化構造";
  return note;
}

function structureLabel(s) {
  if (s.role==="reactant") return "反応前 R";
  if (s.role==="ts") return "認証した TS";
  if (s.role==="product") return /symmetry/.test(s.observed_role || "")?"反応後 P（対称像）":"反応後 P";
  if (s.role==="declared_reactant") return "neutral（宣言した入力）";
  if (s.role==="declared_product") return "shared_proton（宣言した入力）";
  if (s.role==="minimum") return "到達した同じ DFT 極小";
  if (/^monomer/.test(s.role)) return s.atoms.some(a=>a.element==="B")?"BH₃単量体":"NH₃単量体";
  return s.label || s.role;
}

function defaultIndices(structures,kind) {
  const find = pattern => structures.findIndex(s=>pattern.test(s.role || ""));
  if (kind==="tma") return [find(/^declared_reactant$/),find(/^declared_product$/),find(/^minimum$/)].map(i=>Math.max(0,i));
  if (kind==="bh3") return [find(/^monomer_0$/),find(/^monomer_1$/),find(/^product$/)].map(i=>Math.max(0,i));
  const first=find(/reactant|before|input_neutral/);
  const middle=find(/^(ts|saddle)$|transition/);
  const last=find(/product|after|minimum|basin/);
  return [first>=0?first:0,middle>=0?middle:Math.min(1,structures.length-1),last>=0?last:Math.min(2,structures.length-1)];
}

function showMolecules(ex) {
  const structures=validStructures(ex);
  if (!structures.length) return;
  defaultIndices(structures,caseKind(ex)).forEach((value,index)=>{
    // Keep the three WebGL contexts when switching examples.
    if (!$("molecule-"+index)) {
      const card=document.createElement("article"); card.className="molecule-card";
      card.innerHTML=`<div class="molecule-top"><select aria-label="構造${index+1}を選択" id="structure-select-${index}"></select></div><div class="molecule-view" id="molecule-${index}"></div><div class="molecule-bottom"><a id="xyz-link-${index}" download>XYZ</a><div id="structure-info-${index}"></div><div id="distance-${index}"></div></div>`;
      $("molecules").appendChild(card);
    }
    const selector=$("structure-select-"+index);
    selector.innerHTML=structures.map((s,i)=>`<option value="${i}">${esc(structureLabel(s))}</option>`).join("");
    selector.value=String(value);
    try {
      let record=VIEWERS.find(item=>item.index===index);
      if (!record) {
        record={viewer:makeViewer($("molecule-"+index)),selector,structures,index};
        VIEWERS.push(record);
      }
      record.structures=structures;
      const viewer=record.viewer;
      const update=()=>{
        const structure=structures[Number(selector.value)];
        putStructure(viewer,structure,"distance-"+index);
        record.defaultView=viewer.getView();
        $("structure-info-"+index).textContent=structureInfo(structure);
        $("xyz-link-"+index).href=structure.xyz_file ? `data/${structure.xyz_file}` : "#";
      };
      selector.onchange=update; update();
    } catch(error) { $("molecule-"+index).innerHTML=`<div class="error-note">${esc(error.message)} 静的な3D図は下の科学図で確認できます。</div>`; }
  });
}

function showMetrics(ex) {
  const free=ex.free_energy || {};
  const imaginary=ex.decisions?.imag_cm1;
  const values=[
    ["δG_eff",number(free.dG_eff_kcal),"kcal mol⁻¹ · 比較基準は下表"],
    [caseKind(ex)==="oh"?"ΔG_rxn（錯体間）":"ΔG_rxn",number(free.dG_rxn_kcal),"kcal mol⁻¹"],
    ["TSの虚振動",imaginary===null || imaginary===undefined ? "TSなし" : `${number(Math.abs(imaginary),1)}i`,"cm⁻¹ · 電子状態と接続も検証"],
    ["3Dに表示できる構造",String(validStructures(ex).length),"保存座標・出所を保持"]
  ];
  $("metrics").innerHTML=values.map(([label,value,unit])=>`<div><div class="metric-label">${esc(label)}</div><div class="metric-value">${esc(value)}</div><div class="metric-unit">${esc(unit)}</div></div>`).join("");
}

function showEnergy(ex) {
  const free=ex.free_energy;
  if (!free) { $("energy-details").innerHTML="<p>認証した比較点がないため、自由エネルギーの障壁図を作りません。</p>"; return; }
  const points=free.points || [];
  $("energy-details").innerHTML=`<p><strong>比較基準：</strong>${esc(text(free.zero || free.reference))}</p><p><strong>手法：</strong>${esc(text(free.method))}</p><p><strong>条件：</strong>${esc(text(free.conditions))}</p>${table(["保存点","役割","G / Eh","相対G / kcal mol⁻¹"],points.map(p=>[p.label || p.subject,p.role,number(p.G_hartree,8),number(p.relative_kcal,4)]))}<p class="small-note">線の間の自由エネルギーは未計算です。値の欠落はゼロに置き換えません。</p>`;
  if (free.blockers?.length) $("energy-details").innerHTML+=`<p>未認証・不足：${esc(text(free.blockers))}</p>`;
}

function showDecisions(ex,copy) {
  $("decision-flow").innerHTML=copy.flow.map(([name,note],i)=>`${i?'<span class="decision-arrow" aria-hidden="true">→</span>':""}<div class="decision-node"><strong>${esc(name)}</strong>${esc(note)}</div>`).join("");
  $("conclusion").textContent=copy.conclusion; $("limits").textContent=copy.limits;
  const d=ex.decisions || {};
  const facts=[["保存された判定",ex.reaction?.outcome || d.connection_label],["端点接続",d.connection_label],["認証の理由",ex.reaction?.reasons || d.reasons],["虚振動 / cm⁻¹",d.imag_cm1===null || d.imag_cm1===undefined?"該当なし":number(d.imag_cm1,2)],["反応方向のχ",d.chi===null || d.chi===undefined?(d.chi_note || "未記録"):number(d.chi,3)],["blocker",d.blockers?.length?text(d.blockers):"なし（主表示の反応）"]];
  $("decision-facts").innerHTML=facts.map(([label,value])=>`<div class="fact-row"><span>${esc(label)}</span><span>${esc(text(value))}</span></div>`).join("");
  const log=d.log || [];
  $("decision-log").innerHTML=`<ol class="decision-log">${log.map(row=>`<li>${esc(typeof row==="string"?row:JSON.stringify(row))}</li>`).join("")}</ol>`;
  const coverage=ex.coverage;
  $("coverage-section").hidden=!coverage;
  if (coverage) {
    const counts=Array.isArray(coverage.counts)?coverage.counts:Object.values(coverage.counts || {});
    $("coverage").innerHTML=`<p>${esc(text(coverage.grain))}</p><p class="small-note">report時点の保存分類で、DFT入口の選抜結果を含みます。S10はNH₃・HFとTMA・HFのpilot全体。試行済み分類と未試行を集計し、計算失敗はfailure_kindへ別集計します。S10の実NT2試行は3000件、850候補辺から1辺をDFT入口で保持。理由・failureの行は重複する場合があるため、全行を足して成功率にはしません。</p>${table(["探索手順","保存分類の総数","試行済み分類","結果の内訳"],counts.map(row=>[row.mechanism,row.registered_queries,row.executed_queries,JSON.stringify(row.outcomes)]))}<details><summary>保存されたcoverage行</summary>${table(["行"],(coverage.rows || []).map(row=>[JSON.stringify(row)]))}</details>`;
  }
}

function profileStructure(point,ex) {
  if (point.atoms) return {atoms:point.atoms,role:"scan_node",label:point.label || `node ${point.x}`};
  return validStructures(ex).find(s=>s.xyz_file===point.geometry_file || s.role===point.geometry_role);
}

function showProfiles(ex) {
  const profiles=(ex.electronic_profiles || []).filter(p=>p.points?.some(point=>profileStructure(point,ex)));
  $("profile-section").hidden=!profiles.length;
  if (!profiles.length) return;
  $("profile-select").innerHTML=profiles.map((p,i)=>`<option value="${i}">${esc(p.label || p.coordinate_kind || `保存経路 ${i+1}`)}</option>`).join("");
  if (!profileViewer) {
    try { profileViewer=makeViewer($("profile-view")); }
    catch(error) { $("profile-view").innerHTML=`<div class="error-note">${esc(error.message)}</div>`; }
  }
  const pointUpdate=()=>{
    const index=Number($("profile-step").value), point=currentProfile.points[index];
    $("profile-point").textContent=`節点 ${index+1} / ${currentProfile.points.length} · x=${point.x} · ΔE=${number(point.relative_kcal,3)} kcal mol⁻¹`;
    const structure=profileStructure(point,ex);
    if (structure && profileViewer) putStructure(profileViewer,structure,"profile-distance");
    else if (profileViewer) { profileViewer.removeAllModels(); profileViewer.render(); $("profile-distance").textContent="この節点の座標は閲覧用に保存されていません。"; }
  };
  const profileUpdate=()=>{
    currentProfile=profiles[Number($("profile-select").value)];
    $("profile-step").max=String(currentProfile.points.length-1); $("profile-step").value="0";
    $("profile-method").textContent=`${text(currentProfile.method)} / ${text(currentProfile.coordinate_kind)} / zero: ${text(currentProfile.zero)}`;
    $("profile-table").innerHTML=table(["節点","座標","相対E / kcal mol⁻¹"],currentProfile.points.map((p,i)=>[i+1,p.x,number(p.relative_kcal,4)]));
    pointUpdate();
  };
  $("profile-select").onchange=profileUpdate; $("profile-step").oninput=pointUpdate;
  profileUpdate(); if (profileViewer) profileViewer.resize();
}

function showProvenance(ex) {
  const p=ex.provenance || {};
  const runs=(p.runs || []).map(run=>run.run_path || run.run || run.path || run.run_dir || run.root || run);
  $("provenance").innerHTML=`<p><strong>元のrun：</strong><code>${esc(text(runs))}</code></p><p class="small-note">元runを読み取り、図用の小さなデータと構造だけを保存しました。source_filesのSHA-256で出所を照合できます。</p><pre style="white-space:pre-wrap;font-size:11px;overflow-wrap:anywhere">${esc(JSON.stringify(p,null,2))}</pre>`;
}

function selectExample(ex) {
  selectedExample=ex; const kind=caseKind(ex),copy=COPY[kind];
  $("case-kicker").textContent=`EXAMPLE ${String(EXAMPLES.indexOf(ex)+1).padStart(2,"0")} / ${ex.system || ex.id}`;
  $("case-title").textContent=copy.title; $("case-purpose").textContent=copy.purpose;
  $("case-outcome").textContent=copy.outcome; $("case-outcome").className="outcome"+(kind==="tma"?" unresolved":"");
  $("structure-note").textContent=copy.structureNote;
  document.querySelectorAll("#example-nav button").forEach(button=>button.setAttribute("aria-pressed",String(button.dataset.id===ex.id)));
  $("case-figure").src=`figures/${ex.id}.svg`; $("case-figure").alt=`${copy.short}の保存結果：エネルギー、3D構造、判断根拠`;
  $("figure-downloads").innerHTML=`<a href="figures/${esc(ex.id)}.svg" download>SVG</a><a href="figures/${esc(ex.id)}.png" download>PNG 300 dpi</a>`;
  $("figure-caption").textContent=`${copy.short}。図の点と構造は保存済み結果。横軸とゼロ、計算手法は各パネルに表示しています。`;
  showMetrics(ex); showMolecules(ex); showEnergy(ex); showDecisions(ex,copy); showProfiles(ex); showProvenance(ex);
}

$("example-nav").innerHTML=EXAMPLES.map(ex=>`<button type="button" data-id="${esc(ex.id)}" aria-pressed="false">${esc(COPY[caseKind(ex)].short)}</button>`).join("");
$("example-nav").addEventListener("click",event=>{const button=event.target.closest("button[data-id]"); if (button) selectExample(EXAMPLES.find(ex=>ex.id===button.dataset.id));});
$("atom-labels").addEventListener("change",()=>{for (const record of VIEWERS) putStructure(record.viewer,record.structures[Number(record.selector.value)],"distance-"+record.index); if (currentProfile && profileViewer) $("profile-step").oninput();});
$("reset-view").addEventListener("click",()=>{for (const {viewer,defaultView} of VIEWERS) {viewer.setView(defaultView);viewer.render();}});
$("sync-view").addEventListener("click",()=>{if (!VIEWERS.length) return; const source=VIEWERS[0].viewer.getView(); for (const {viewer} of VIEWERS.slice(1)) {const view=viewer.getView();view.splice(4,4,...source.slice(4,8));viewer.setView(view);viewer.render();}});
window.addEventListener("resize",()=>{for (const {viewer} of VIEWERS) viewer.resize(); if (profileViewer && !$("profile-section").hidden) profileViewer.resize();});

const references=Array.isArray(PAYLOAD.references)?PAYLOAD.references:(PAYLOAD.references.references || []);
$("references").innerHTML=references.map(ref=>`<li id="ref-${esc(ref.id)}"><a href="${esc(ref.url)}" target="_blank" rel="noopener">${esc(ref.citation || ref.title)}</a><div class="ref-use">${esc(ref.role || "")}</div>${(ref.related_papers || []).map(r=>`<div class="ref-use">関連：<a href="${esc(r.url || `https://doi.org/${r.doi}`)}" target="_blank" rel="noopener">${esc(r.citation || r.title || r.doi)}</a></div>`).join("")}</li>`).join("");
$("library-note").textContent=`表示：3Dmol.js ${PAYLOAD.library.version} (${PAYLOAD.library.license}) を同梱。静的図：Matplotlib。データ読込と3D表示にネット接続は不要です。`;
selectExample(EXAMPLES.find(ex=>caseKind(ex)==="inversion") || EXAMPLES[0]);
