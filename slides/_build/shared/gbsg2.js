/* GBSG2 breast cancer trial (lifelines): the "population" used by the Day 1 simulations. Inlined by build.py.
   Needs deck.js first (mean, sdP, skew). */
const DATA = __DATA__;
const VARS = {
  tsize:{label:'Tumour size', unit:'mm', bw:5, dom:[0,120]},
  age:{label:'Age', unit:'years', bw:2, dom:[20,82]},
  pnodes:{label:'Positive lymph nodes', unit:'nodes', bw:1, dom:[0,52]}
};
const POP = {}; for(const k in VARS){const a=DATA[k];POP[k]={mu:mean(a),sd:sdP(a),med:d3.median(a),skew:skew(a),min:d3.min(a),max:d3.max(a)}}
