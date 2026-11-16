/* GIMM interactive deck: shared engine and helpers. Inlined into every deck by build.py.
   A deck defines  const hooks = { "<slide id>": {init(), enter(), leave()} }  and calls  startDeck(hooks). */
const C = {pc1:'#240046', pc2:'#9933ff', sc1:'#ff0066', sc3:'#ff99cc', pop:'#a799b5', popL:'#d3cce0', ink2:'#5b4b6e', line:'#e6e1ec'};
/* ---------- maths ---------- */
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296}}
let R = mulberry32((Date.now()%1e9)|0);
const mean = a => {let s=0;for(const v of a)s+=v;return s/a.length};
const sdP = a => {const m=mean(a);let s=0;for(const v of a)s+=(v-m)*(v-m);return Math.sqrt(s/a.length)};
const sdS = a => {const m=mean(a);let s=0;for(const v of a)s+=(v-m)*(v-m);return Math.sqrt(s/(a.length-1))};
const skew = a => {const m=mean(a),s=sdP(a);let t=0;for(const v of a)t+=((v-m)/s)**3;return t/a.length};
function sampleIdx(N,n,r=R){const o=new Array(n);for(let i=0;i<n;i++)o[i]=Math.floor(r()*N);return o}
function studyMeans(pop,n,k,r=R,cum=null){
  const N=pop.length,out=new Float64Array(k);
  for(let j=0;j<k;j++){let s=0;for(let i=0;i<n;i++){let idx;
    if(cum){const u=r()*cum[N-1];let lo=0,hi=N-1;while(lo<hi){const mid=(lo+hi)>>1;if(cum[mid]<u)lo=mid+1;else hi=mid}idx=lo}else idx=Math.floor(r()*N);
    s+=pop[idx]}out[j]=s/n}
  return out;
}
const fmt1 = d3.format('.1f'), fmt2=d3.format('.2f');

/* ---------- deck chrome (progress bar, nav, tooltip, notes) ---------- */
(function(){
  const stage=document.getElementById('stage');
  stage.insertAdjacentHTML('beforeend','<div id="progress"></div><div id="counter"></div><div id="nav"><button id="prevBtn" aria-label="Previous">‹</button><button id="nextBtn" aria-label="Next">›</button></div>');
  document.body.insertAdjacentHTML('beforeend','<div id="tip"></div><div id="notes"><h4>Speaker notes</h4><div id="notesBody"></div></div><div id="help">← → navigate · N notes · F fullscreen</div>');
})();

/* ---------- tooltip ---------- */
const tip=document.getElementById('tip');
function showTip(ev,html){tip.innerHTML=html;tip.style.opacity=1;let x=ev.clientX+14,y=ev.clientY-12;if(x>innerWidth-220)x=ev.clientX-200;tip.style.left=x+'px';tip.style.top=y+'px'}
function hideTip(){tip.style.opacity=0}


/* ---------- small helpers ---------- */
function axisB(g,x,ticks,fmt){g.call(d3.axisBottom(x).ticks(ticks).tickSize(5).tickFormat(fmt||null)).attr('class','axis')}
function axisL(g,y,ticks,w){g.call(d3.axisLeft(y).ticks(ticks).tickSize(-w).tickPadding(8)).attr('class','axis grid').call(s=>s.selectAll('.tick line').attr('stroke',C.line)).call(s=>s.select('.domain').remove())}
function segButtons(container,items,init,onPick){
  const el=document.getElementById(container);const btns={};
  items.forEach(([val,label])=>{const b=document.createElement('button');b.className='btn'+(val===init?' on':'');b.textContent=label;b.onclick=()=>{Object.values(btns).forEach(x=>x.classList.remove('on'));b.classList.add('on');onPick(val)};el.appendChild(b);btns[val]=b});
}
function barPath(x0,x1,y0,yTop,r){ // bar with rounded top corners anchored at baseline y0
  const w=x1-x0,h=y0-yTop; if(h<=0.5||w<=0)return `M${x0},${y0}h${Math.max(w,0)}`;
  r=Math.min(r,w/2,h);return `M${x0},${y0}V${yTop+r}Q${x0},${yTop} ${x0+r},${yTop}H${x1-r}Q${x1},${yTop} ${x1},${yTop+r}V${y0}Z`;
}


/* ---------- deck engine ---------- */
function startDeck(hooks){
  const slides=[...document.querySelectorAll('.slide')];
  const stage=document.getElementById('stage');
  let cur=-1;const inited=new Set();
  function fit(){const s=Math.min(innerWidth/1280,innerHeight/720);stage.style.transform=`scale(${s})`}
  addEventListener('resize',fit);fit();
  const notesEl=document.getElementById('notes'),notesBody=document.getElementById('notesBody');
  function go(i,push=true){
    i=Math.max(0,Math.min(slides.length-1,i));if(i===cur)return;
    if(cur>=0){slides[cur].classList.remove('active');const h=hooks[slides[cur].id];h&&h.leave&&h.leave()}
    cur=i;const s=slides[i];s.classList.add('active');hideTip();
    const h=hooks[s.id];
    if(h){if(!inited.has(s.id)){inited.add(s.id);try{h.init()}catch(e){console.error(e)}}else if(h.enter)h.enter()}
    document.getElementById('progress').style.width=((i+1)/slides.length*100)+'%';
    document.getElementById('counter').textContent=`${i+1} / ${slides.length}`;
    const n=s.querySelector('aside.notes');notesBody.textContent=n?n.textContent.trim():'';
    if(push){try{history.replaceState(null,'','#'+s.id)}catch(e){location.hash=s.id}}
  }
  function fromHash(){const id=location.hash.slice(1);const i=slides.findIndex(s=>s.id===id);go(i>=0?i:0,false)}
  addEventListener('hashchange',fromHash);
  document.getElementById('prevBtn').onclick=()=>go(cur-1);
  document.getElementById('nextBtn').onclick=()=>go(cur+1);
  addEventListener('keydown',e=>{
    if(e.target.tagName==='INPUT')return;
    if(['ArrowRight','PageDown'].includes(e.key)||(e.key===' '&&e.target.tagName!=='BUTTON')){e.preventDefault();go(cur+1)}
    else if(['ArrowLeft','PageUp'].includes(e.key)){e.preventDefault();go(cur-1)}
    else if(e.key==='Home')go(0);else if(e.key==='End')go(slides.length-1);
    else if(e.key==='n'||e.key==='N')notesEl.classList.toggle('show');
    else if(e.key==='f'||e.key==='F'){try{document.fullscreenElement?document.exitFullscreen():document.documentElement.requestFullscreen()}catch(err){}}
  });
  fromHash();
}
