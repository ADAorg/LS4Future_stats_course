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
// one call = k studies of n patients: each study's mean and sample SD (ddof = 1)
function studyStats(pop,n,k,r=R){const N=pop.length,m=new Float64Array(k),s=new Float64Array(k);
  for(let j=0;j<k;j++){let a=0,b=0;for(let i=0;i<n;i++){const v=pop[Math.floor(r()*N)];a+=v;b+=v*v}const mu=a/n;m[j]=mu;s[j]=Math.sqrt(Math.max(0,(b-n*mu*mu)/(n-1)))}
  return{m,s};
}
/* ---------- distributions (checked against scipy.stats to 5 decimals) ---------- */
function lgamma(x){const c=[0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61503916999185,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
  if(x<0.5)return Math.log(Math.PI/Math.sin(Math.PI*x))-lgamma(1-x);x-=1;let a=c[0];const t=x+7.5;for(let i=1;i<9;i++)a+=c[i]/(x+i);return 0.5*Math.log(2*Math.PI)+(x+0.5)*Math.log(t)-t+Math.log(a)}
function betacf(a,b,x){const E=3e-14,F=1e-300;let qab=a+b,qap=a+1,qam=a-1,c=1,d=1-qab*x/qap;if(Math.abs(d)<F)d=F;d=1/d;let h=d;
  for(let m=1;m<=300;m++){const m2=2*m;let aa=m*(b-m)*x/((qam+m2)*(a+m2));d=1+aa*d;if(Math.abs(d)<F)d=F;c=1+aa/c;if(Math.abs(c)<F)c=F;d=1/d;h*=d*c;
    aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));d=1+aa*d;if(Math.abs(d)<F)d=F;c=1+aa/c;if(Math.abs(c)<F)c=F;d=1/d;const del=d*c;h*=del;if(Math.abs(del-1)<E)break}return h}
function ibeta(x,a,b){if(x<=0)return 0;if(x>=1)return 1;const bt=Math.exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*Math.log(x)+b*Math.log(1-x));return x<(a+1)/(a+b+2)?bt*betacf(a,b,x)/a:1-bt*betacf(b,a,1-x)/b}
function tCDF(t,df){const p=0.5*ibeta(df/(df+t*t),df/2,0.5);return t>0?1-p:p}
function tQuantile(p,df){let lo=-1e3,hi=1e3;for(let i=0;i<200;i++){const m=(lo+hi)/2;if(tCDF(m,df)<p)lo=m;else hi=m}return(lo+hi)/2}
const normPDF=(x,mu=0,sd=1)=>Math.exp(-0.5*((x-mu)/sd)**2)/(sd*Math.sqrt(2*Math.PI));
function binomPMF(k,n,p){return Math.exp(lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)+k*Math.log(p)+(n-k)*Math.log(1-p))}
function binomQuantile(q,n,p){let c=0;for(let k=0;k<=n;k++){c+=binomPMF(k,n,p);if(c>=q)return k}return n}

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
