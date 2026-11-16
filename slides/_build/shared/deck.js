/* GIMM interactive deck: shared engine and helpers. Inlined into every deck by build.py.
   A deck defines  const hooks = { "<slide id>": {init(), enter(), leave()} }  and calls  startDeck(hooks). */
const C = {pc1:'#240046', pc2:'#9933ff', sc1:'#ff0066', sc3:'#ff99cc', pop:'#a799b5', popL:'#d3cce0', ink2:'#5b4b6e', muted:'#8a7d99', line:'#e6e1ec'};
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
function normCDF(x){const t=1/(1+0.3275911*Math.abs(x)/Math.SQRT2),y=1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-0.284496736)*t+0.254829592)*t*Math.exp(-x*x/2);return x>=0?(1+y)/2:(1-y)/2}  // |error| < 1e-7
function normQuantile(p){let lo=-40,hi=40;for(let i=0;i<100;i++){const m=(lo+hi)/2;if(normCDF(m)<p)lo=m;else hi=m}return(lo+hi)/2}
const normPDF=(x,mu=0,sd=1)=>Math.exp(-0.5*((x-mu)/sd)**2)/(sd*Math.sqrt(2*Math.PI));
function binomPMF(k,n,p){return Math.exp(lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)+k*Math.log(p)+(n-k)*Math.log(1-p))}
function binomQuantile(q,n,p){let c=0;for(let k=0;k<=n;k++){c+=binomPMF(k,n,p);if(c>=q)return k}return n}
/* ---------- tests (checked against scipy.stats) ---------- */
const varS = a => {const m=mean(a);let s=0;for(const v of a)s+=(v-m)*(v-m);return s/(a.length-1)};
function tTwoSided(t,df){return ibeta(df/(df+t*t),df/2,0.5)}   // exact two-sided tail, no 1-x precision loss
function welch(a,b){const ma=mean(a),mb=mean(b),va=varS(a)/a.length,vb=varS(b)/b.length,se=Math.sqrt(va+vb);
  const df=(va+vb)**2/(va*va/(a.length-1)+vb*vb/(b.length-1)),t=(ma-mb)/se;return{diff:ma-mb,se,df,t,p:tTwoSided(Math.abs(t),df)}}
function gammaQ(a,x){ // regularized upper incomplete gamma
  if(x<=0)return 1;const gln=lgamma(a);
  if(x<a+1){let ap=a,sum=1/a,del=sum;for(let n=0;n<500;n++){ap++;del*=x/ap;sum+=del;if(Math.abs(del)<Math.abs(sum)*1e-15)break}return 1-sum*Math.exp(-x+a*Math.log(x)-gln)}
  let b=x+1-a,c=1e300,d=1/b,h=d;for(let i=1;i<500;i++){const an=-i*(i-a);b+=2;d=an*d+b;if(Math.abs(d)<1e-300)d=1e-300;c=b+an/c;if(Math.abs(c)<1e-300)c=1e-300;d=1/d;const del=d*c;h*=del;if(Math.abs(del-1)<1e-15)break}
  return Math.exp(-x+a*Math.log(x)-gln)*h}
const chi2SF=(x,k)=>gammaQ(k/2,x/2);
function chi2Table(g,c){ // g: 0/1 group codes, c: category codes; Yates correction for 2x2 (as scipy.stats.chi2_contingency)
  const cats=[...new Set(c)].sort((a,b)=>a-b),gs=[0,1],O=gs.map(gg=>cats.map(cc=>0));
  for(let i=0;i<g.length;i++)O[g[i]][cats.indexOf(c[i])]++;
  const N=g.length,rs=O.map(r=>d3.sum(r)),cs=cats.map((_,j)=>O[0][j]+O[1][j]),df=(cats.length-1);let x=0;
  for(let i=0;i<2;i++)for(let j=0;j<cats.length;j++){const E=rs[i]*cs[j]/N;let d=Math.abs(O[i][j]-E);if(df===1)d=Math.max(0,d-0.5);x+=d*d/E}
  return{stat:x,df,p:chi2SF(x,df),O,cats,rs}}
function randn(r=R){let u=0;while(u===0)u=r();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*r())}
// k label shuffles: difference in means between a random group of size nA and the rest
function permDiffs(vals,nA,k,r=R){const N=vals.length,idx=Int32Array.from({length:N},(_,i)=>i);let tot=0;for(const v of vals)tot+=v;const out=new Float64Array(k);
  for(let j=0;j<k;j++){let s=0;for(let i=0;i<nA;i++){const q=i+Math.floor(r()*(N-i)),t=idx[i];idx[i]=idx[q];idx[q]=t;s+=vals[idx[i]]}out[j]=s/nA-(tot-s)/(N-nA)}return out}
function shuffled(a,r=R){const b=a.slice();for(let i=b.length-1;i>0;i--){const j=Math.floor(r()*(i+1));[b[i],b[j]]=[b[j],b[i]]}return b}
function fmtP(p){if(p>=0.001)return d3.format('.2~g')(p).replace(/^0\./,'0.');const e=Math.floor(Math.log10(p)),m=p/10**e;
  const sup=String(e).replace('-','⁻').replace(/\d/g,d=>'⁰¹²³⁴⁵⁶⁷⁸⁹'[d]);return `${m.toFixed(1)} × 10${sup}`}
/* ---------- more tests: paired, rank-based, exact, survival (checked against scipy / lifelines) ---------- */
function rankAvg(v){const o=v.map((x,i)=>[x,i]).sort((a,b)=>a[0]-b[0]),r=new Array(v.length);let ties=0;
  for(let i=0;i<o.length;){let j=i;while(j+1<o.length&&o[j+1][0]===o[i][0])j++;const rk=(i+j)/2+1;for(let k=i;k<=j;k++)r[o[k][1]]=rk;const t=j-i+1;ties+=t*t*t-t;i=j+1}
  return{r,ties}}
function pairedT(x,y){const d=x.map((v,i)=>v-y[i]),n=d.length,md=mean(d),se=Math.sqrt(varS(d)/n),t=md/se;return{diff:md,se,df:n-1,t,p:tTwoSided(Math.abs(t),n-1)}}
function wilcoxonSR(x,y){ // exact two-sided (no ties, no zeros, n <= 25), like scipy's default for small n
  const d=x.map((v,i)=>v-y[i]).filter(v=>v!==0),n=d.length,{r}=rankAvg(d.map(Math.abs));let wp=0;d.forEach((v,i)=>{if(v>0)wp+=r[i]});
  const M=n*(n+1)/2,cnt=new Float64Array(M+1);cnt[0]=1;for(let k=1;k<=n;k++)for(let s=M;s>=k;s--)cnt[s]+=cnt[s-k];
  const tot=2**n;let lo=0,hi=0;for(let s=0;s<=M;s++){if(s<=wp)lo+=cnt[s];if(s>=wp)hi+=cnt[s]}
  return{W:Math.min(wp,M-wp),p:Math.min(1,2*Math.min(lo,hi)/tot)}}
function mannWhitney(x,y){ // exact when both n <= 8 and no ties, else normal approximation with tie and continuity correction (scipy 'auto')
  const n1=x.length,n2=y.length,all=x.concat(y),{r,ties}=rankAvg(all);let R1=0;for(let i=0;i<n1;i++)R1+=r[i];
  const U1=R1-n1*(n1+1)/2,U=Math.max(U1,n1*n2-U1);
  if(n1<=8&&n2<=8&&ties===0){ // count arrangements with U >= observed
    const f=(()=>{const memo=new Map();const g=(a,b,u)=>{if(u<0)return 0;if(a===0||b===0)return u===0?1:0;const k=a+','+b+','+u;if(memo.has(k))return memo.get(k);const v=g(a-1,b,u-b)+g(a,b-1,u);memo.set(k,v);return v};return g})();
    let tot=0,ge=0;for(let u=0;u<=n1*n2;u++){const c=f(n1,n2,u);tot+=c;if(u>=U)ge+=c}return{U:U1,p:Math.min(1,2*ge/tot)}}
  const mu=n1*n2/2,N=n1+n2,s=Math.sqrt(n1*n2/12*((N+1)-ties/(N*(N-1))));const z=(U-mu-0.5)/s;return{U:U1,p:Math.min(1,2*(1-normCDF(z)))}}
function fisher2x2(t){ // t = [[a,b],[c,d]], two-sided
  const [a,b]=t[0],[c,d]=t[1],r1=a+b,c1=a+c,N=a+b+c+d;const lc=(n,k)=>lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1);
  const pr=k=>Math.exp(lc(c1,k)+lc(N-c1,r1-k)-lc(N,r1));const po=pr(a);let p=0;
  for(let k=Math.max(0,r1-(N-c1));k<=Math.min(r1,c1);k++){const q=pr(k);if(q<=po*(1+1e-7))p+=q}return Math.min(1,p)}
function chi2Plain(t){const R=t.map(r=>d3.sum(r)),Cc=t[0].map((_,j)=>d3.sum(t,r=>r[j])),N=d3.sum(R);let x=0;const E=t.map((r,i)=>r.map((_,j)=>R[i]*Cc[j]/N));
  t.forEach((r,i)=>r.forEach((o,j)=>x+=(o-E[i][j])**2/E[i][j]));const df=(t.length-1)*(t[0].length-1);return{stat:x,df,p:chi2SF(x,df),E}}
function kaplanMeier(time,event){const idx=d3.range(time.length).sort((a,b)=>time[a]-time[b]);let atRisk=time.length,S=1;const pts=[[0,1]];
  for(let i=0;i<idx.length;){const t=time[idx[i]];let d=0,c=0;while(i<idx.length&&time[idx[i]]===t){if(event[idx[i]])d++;else c++;i++}
    if(d>0){S*=1-d/atRisk;pts.push([t,S])}atRisk-=d+c}return pts}
function logRank(t1,e1,t2,e2){const T=[...new Set([...t1.filter((_,i)=>e1[i]),...t2.filter((_,i)=>e2[i])])].sort((a,b)=>a-b);let O=0,E=0,V=0;
  for(const t of T){let n1=0,n2=0,d1=0,d2=0;t1.forEach((v,i)=>{if(v>=t){n1++;if(v===t&&e1[i])d1++}});t2.forEach((v,i)=>{if(v>=t){n2++;if(v===t&&e2[i])d2++}});
    const n=n1+n2,d=d1+d2;O+=d1;E+=d*n1/n;if(n>1)V+=d*(n1/n)*(n2/n)*(n-d)/(n-1)}
  const stat=(O-E)**2/V;return{stat,p:chi2SF(stat,1),O,E}}
/* ---------- linear model (OLS), checked against statsmodels ---------- */
function matInv(A){const n=A.length,M=A.map((r,i)=>[...r,...A.map((_,j)=>i===j?1:0)]);
  for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(M[r][c])>Math.abs(M[p][c]))p=r;[M[c],M[p]]=[M[p],M[c]];
    const v=M[c][c];for(let j=0;j<2*n;j++)M[c][j]/=v;for(let r=0;r<n;r++)if(r!==c){const f=M[r][c];if(f)for(let j=0;j<2*n;j++)M[r][j]-=f*M[c][j]}}
  return M.map(r=>r.slice(n))}
const fSF=(F,d1,d2)=>F<=0?1:ibeta(d2/(d2+d1*F),d2/2,d1/2);
// X: array of rows (include a leading 1 for the intercept); returns coefficients with SE, t, p, 95% CI, and the overall F-test
function ols(X,y){const n=y.length,k=X[0].length,XtX=Array.from({length:k},()=>new Array(k).fill(0)),Xty=new Array(k).fill(0);
  for(let i=0;i<n;i++){const r=X[i];for(let a=0;a<k;a++){Xty[a]+=r[a]*y[i];for(let b=0;b<k;b++)XtX[a][b]+=r[a]*r[b]}}
  const V=matInv(XtX),b=V.map(row=>row.reduce((s,v,j)=>s+v*Xty[j],0));let rss=0;const my=mean(y);let tss=0;const fitted=new Array(n);
  for(let i=0;i<n;i++){const f=X[i].reduce((s,v,j)=>s+v*b[j],0);fitted[i]=f;rss+=(y[i]-f)**2;tss+=(y[i]-my)**2}
  const df=n-k,s2=rss/df,se=V.map((r,i)=>Math.sqrt(r[i]*s2)),tc=tQuantile(.975,df);
  const coef=b.map((v,i)=>({est:v,se:se[i],t:v/se[i],p:tTwoSided(Math.abs(v/se[i]),df),lo:v-tc*se[i],hi:v+tc*se[i]}));
  const F=((tss-rss)/(k-1))/s2;return{coef,df,rss,sigma:Math.sqrt(s2),r2:1-rss/tss,F,Fp:k>1?fSF(F,k-1,df):NaN,fitted,resid:y.map((v,i)=>v-fitted[i])}}
// design matrix helpers
function dummies(vals,levels,ref){return vals.map(v=>levels.filter(l=>l!==ref).map(l=>v===l?1:0))}
/* ---------- regression diagnostics (checked against statsmodels) ---------- */
// leverage, internally studentized residuals, Cook's distance, HC3 standard errors
function olsDiag(X,y){const f=ols(X,y),n=y.length,k=X[0].length,XtX=Array.from({length:k},()=>new Array(k).fill(0));
  X.forEach(r=>{for(let a=0;a<k;a++)for(let b=0;b<k;b++)XtX[a][b]+=r[a]*r[b]});const V=matInv(XtX);
  const h=X.map(r=>{let s=0;for(let a=0;a<k;a++)for(let b=0;b<k;b++)s+=r[a]*V[a][b]*r[b];return s});
  const s2=f.rss/f.df,std=f.resid.map((e,i)=>e/Math.sqrt(s2*(1-h[i]))),cook=std.map((r,i)=>r*r/k*h[i]/(1-h[i]));
  const M=Array.from({length:k},()=>new Array(k).fill(0));X.forEach((r,i)=>{const w=(f.resid[i]/(1-h[i]))**2;for(let a=0;a<k;a++)for(let b=0;b<k;b++)M[a][b]+=w*r[a]*r[b]});
  const VM=V.map(row=>row.map((_,j)=>row.reduce((s,v,l)=>s+v*M[l][j],0))),HC=VM.map(row=>row.map((_,j)=>row.reduce((s,v,l)=>s+v*V[l][j],0)));
  f.h=h;f.std=std;f.cook=cook;f.seHC3=HC.map((r,i)=>Math.sqrt(r[i]));return f}
// LOWESS smoother (tricube weights, local linear, frac of points, robustness iterations as statsmodels' default it=3)
function lowess(xs,ys,frac=2/3,it=3){const n=xs.length,o=d3.range(n).sort((a,b)=>xs[a]-xs[b]),x=o.map(i=>xs[i]),y=o.map(i=>ys[i]),r=Math.max(2,Math.ceil(frac*n));let rw=new Array(n).fill(1),fit=new Array(n);
  for(let t=0;t<=it;t++){for(let i=0;i<n;i++){const d=x.map(v=>Math.abs(v-x[i])),hmax=[...d].sort((a,b)=>a-b)[r-1]||1e-12;let sw=0,sx=0,sy=0,sxx=0,sxy=0;
      for(let j=0;j<n;j++){const u=d[j]/hmax,w=u<1?(1-u**3)**3*rw[j]:0;sw+=w;sx+=w*x[j];sy+=w*y[j];sxx+=w*x[j]*x[j];sxy+=w*x[j]*y[j]}
      const den=sw*sxx-sx*sx;fit[i]=Math.abs(den)<1e-12?sy/sw:((sxx*sy-sx*sxy)+(sw*sxy-sx*sy)*x[i])/den}
    if(t===it)break;const res=y.map((v,i)=>Math.abs(v-fit[i])),med=[...res].sort((a,b)=>a-b)[Math.floor(n/2)]||1e-12;rw=res.map(e=>{const u=e/(6*med);return u<1?(1-u*u)**2:0})}
  return x.map((v,i)=>[v,fit[i]])}
/* ---------- logistic regression (IRLS), checked against statsmodels ---------- */
function logitFit(X,y,iter=30){const n=y.length,k=X[0].length;let b=new Array(k).fill(0),V;
  for(let t=0;t<iter;t++){const A=Array.from({length:k},()=>new Array(k).fill(0)),g=new Array(k).fill(0);
    for(let i=0;i<n;i++){const r=X[i],eta=r.reduce((s,v,j)=>s+v*b[j],0),p=1/(1+Math.exp(-eta)),w=Math.max(p*(1-p),1e-10);
      for(let a=0;a<k;a++){g[a]+=r[a]*(y[i]-p);for(let c=0;c<k;c++)A[a][c]+=w*r[a]*r[c]}}
    V=matInv(A);const step=V.map(row=>row.reduce((s,v,j)=>s+v*g[j],0));b=b.map((v,j)=>v+step[j]);if(Math.max(...step.map(Math.abs))<1e-10)break}
  const coef=b.map((v,i)=>{const se=Math.sqrt(V[i][i]),z=v/se;return{est:v,se,z,p:chi2SF(z*z,1),or:Math.exp(v),lo:Math.exp(v-1.959964*se),hi:Math.exp(v+1.959964*se)}});
  const pred=r=>1/(1+Math.exp(-r.reduce((s,v,j)=>s+v*b[j],0)));return{coef,b,pred}}
/* ---------- step widget: numbered dots + Back / Next inside a container ---------- */
function stepper(id,n,onStep){const el=document.getElementById(id);let cur=0;const dots=[];
  const prev=document.createElement('button');prev.className='btn';prev.textContent='‹ Back';const next=document.createElement('button');next.className='btn primary';next.textContent='Next step ›';
  for(let i=0;i<n;i++){const d=document.createElement('div');d.className='dot';d.textContent=i+1;d.onclick=()=>go(i);el.appendChild(d);dots.push(d)}
  el.appendChild(prev);el.appendChild(next);
  function go(i){cur=Math.max(0,Math.min(n-1,i));dots.forEach((d,k)=>{d.className='dot'+(k===cur?' on':(k<cur?' done':''))});prev.disabled=cur===0;next.disabled=cur===n-1;onStep(cur)}
  prev.onclick=()=>go(cur-1);next.onclick=()=>go(cur+1);return{go}}
/* ---------- concept-check cards: Q = [[isCorrect, statement, explanation], ...] ---------- */
function qcards(id,Q){const box=document.getElementById(id);
  Q.forEach(([ok,q,a],i)=>{const d=document.createElement('div');d.className='qcard '+(ok?'ok':'no');
    d.innerHTML=`<div class="tag">${i+1} · <span class="res" style="visibility:hidden">${ok?'✓ correct':'✗ wrong'}</span></div><div class="q">${q}</div><div class="a">${a}</div>`;
    const res=d.querySelector('.res');d.onclick=()=>{d.classList.toggle('shown');res.style.visibility=d.classList.contains('shown')?'visible':'hidden'};box.appendChild(d)});}

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
