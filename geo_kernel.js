/* ===== 立体几何 JS 计算核心（移植自 edulab geometry_kernel.py，输出 edulab 数据 schema） =====
   坐标为数学坐标(x,y,z)，z 轴竖直；模板 threeFromMath 把 (x,y,z)->three(x,z,y)。 */
(function(){
  const V={
    sub:(a,b)=>[a[0]-b[0],a[1]-b[1],a[2]-b[2]],
    add:(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]],
    scale:(a,k)=>[a[0]*k,a[1]*k,a[2]*k],
    dot:(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],
    cross:(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],
    len:a=>Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]),
    mid:(a,b)=>[(a[0]+b[0])/2,(a[1]+b[1])/2,(a[2]+b[2])/2]
  };
  function simpRad(N){let k=1,m=Math.round(N);if(m<0)m=0;for(let d=2;d*d<=m;d++){while(m%(d*d)===0){m/=d*d;k*=d;}}return{k:k,m:m};}
  function lenLatex(x){const sq=x*x,N=Math.round(sq);
    if(Math.abs(sq-N)<1e-6){const r=simpRad(N);if(r.m===1)return''+r.k;return (r.k>1?r.k:'')+'\\sqrt{'+r.m+'}';}
    return (+x.toFixed(3)).toString();}
  function fracTrigLatex(x){const ax=Math.abs(x),sgn=x<0?'-':'';
    for(let b=1;b<=12;b++)for(let a=0;a<=b*b;a++){if(Math.abs(ax-Math.sqrt(a)/b)<1e-6){
      if(a===0)return'0';const r=simpRad(a);const num=(r.m===1)?(''+r.k):((r.k>1?r.k:'')+'\\sqrt{'+r.m+'}');
      if(b===1)return sgn+num;return sgn+'\\frac{'+num+'}{'+b+'}';}}
    return (+x.toFixed(4)).toString();}
  function numLatex(x){const r=Math.round(x);return Math.abs(x-r)<1e-9?(''+r):(+x.toFixed(3)).toString();}
  function deg(x){return (+x.toFixed(1)).toString()+'^\\circ';}
  function coordTex(name,p){return name+'('+p.map(numLatex).join(',')+')';}
  function vecTex(name,v){return '\\vec{'+name+'}=('+v.map(numLatex).join(',')+')';}
  function threeOf(p){return [p[0],p[2],p[1]];}                 // math -> three 坐标
  function centroid(points){const ns=Object.keys(points);const s=[0,0,0];ns.forEach(n=>{const p=points[n];s[0]+=p[0];s[1]+=p[1];s[2]+=p[2];});return [s[0]/ns.length,s[1]/ns.length,s[2]/ns.length];}
  function maxSpan(points){let mx=1;Object.values(points).forEach(p=>{mx=Math.max(mx,Math.abs(p[0]),Math.abs(p[1]),Math.abs(p[2]));});return mx;}

  /* ---- 几何体（数学坐标，z 竖直） ---- */
  function cube(a){a=a||2;return{
    points:{A:[0,0,0],B:[a,0,0],C:[a,a,0],D:[0,a,0],A1:[0,0,a],B1:[a,0,a],C1:[a,a,a],D1:[0,a,a]},
    edges:[['A','B'],['B','C'],['C','D'],['D','A'],['A1','B1'],['B1','C1'],['C1','D1'],['D1','A1'],['A','A1'],['B','B1'],['C','C1'],['D','D1']]};}
  function cuboid(lx,ly,lz){lx=lx||3;ly=ly||2;lz=lz||2;return{
    points:{A:[0,0,0],B:[lx,0,0],C:[lx,ly,0],D:[0,ly,0],A1:[0,0,lz],B1:[lx,0,lz],C1:[lx,ly,lz],D1:[0,ly,lz]},
    edges:[['A','B'],['B','C'],['C','D'],['D','A'],['A1','B1'],['B1','C1'],['C1','D1'],['D1','A1'],['A','A1'],['B','B1'],['C','C1'],['D','D1']]};}
  function pyramid4(a,h){a=a||2;h=h||3;return{
    points:{A:[0,0,0],B:[a,0,0],C:[a,a,0],D:[0,a,0],P:[a/2,a/2,h]},
    edges:[['A','B'],['B','C'],['C','D'],['D','A'],['A','P'],['B','P'],['C','P'],['D','P']]};}
  function tetra(a){a=a||2;const h=a*Math.sqrt(2/3);return{
    points:{A:[0,0,0],B:[a,0,0],C:[a/2,a*Math.sqrt(3)/2,0],D:[a/2,a*Math.sqrt(3)/6,h]},
    edges:[['A','B'],['B','C'],['C','A'],['A','D'],['B','D'],['C','D']]};}
  function prism3(a,h){a=a||2;h=h||3;const s=a*Math.sqrt(3)/2;return{
    points:{A:[0,0,0],B:[a,0,0],C:[a/2,s,0],A1:[0,0,h],B1:[a,0,h],C1:[a/2,s,h]},
    edges:[['A','B'],['B','C'],['C','A'],['A1','B1'],['B1','C1'],['C1','A1'],['A','A1'],['B','B1'],['C','C1']]};}

  function tokPts(run){return run.match(/[A-Z]\d?/g)||[];}
  function normSub(t){return (t||'').replace(/[₀₁₂₃₄₅₆₇₈₉]/g,function(c){return '₀₁₂₃₄₅₆₇₈₉'.indexOf(c);});}
  function num(t,re,d){const m=t.match(re);return m?+m[1]:d;}

  function detectBody(t){
    if(/正方体|正四棱柱/.test(t)){const a=num(t,/棱长[为是=]?\s*(\d+\.?\d*)/)||num(t,/边长[为是=]?\s*(\d+\.?\d*)/)||2;return{key:'cube',geo:cube(a),title:'正方体 ABCD-A₁B₁C₁D₁（棱长 '+a+'）',a:a};}
    if(/长方体/.test(t)){const lx=num(t,/长[为是=]?\s*(\d+\.?\d*)/,3),ly=num(t,/宽[为是=]?\s*(\d+\.?\d*)/,2),lz=num(t,/高[为是=]?\s*(\d+\.?\d*)/,2);return{key:'cuboid',geo:cuboid(lx,ly,lz),title:'长方体 ABCD-A₁B₁C₁D₁（'+lx+'×'+ly+'×'+lz+'）'};}
    if(/正四棱锥|四棱锥/.test(t)){const a=num(t,/底面边长[为是=]?\s*(\d+\.?\d*)/)||num(t,/底边[为是=]?\s*(\d+\.?\d*)/)||2,h=num(t,/高[为是=]?\s*(\d+\.?\d*)/,3);return{key:'pyramid4',geo:pyramid4(a,h),title:'正四棱锥 P-ABCD（底边 '+a+'，高 '+h+'）'};}
    if(/正四面体|正三棱锥/.test(t)){const a=num(t,/棱长[为是=]?\s*(\d+\.?\d*)/)||num(t,/边长[为是=]?\s*(\d+\.?\d*)/)||2;return{key:'tetra',geo:tetra(a),title:'正四面体 ABCD（棱长 '+a+'）'};}
    if(/正三棱柱|三棱柱/.test(t)){const a=num(t,/底面边长[为是=]?\s*(\d+\.?\d*)/)||num(t,/边长[为是=]?\s*(\d+\.?\d*)/)||2,h=num(t,/高[为是=]?\s*(\d+\.?\d*)/)||num(t,/侧棱[为是=]?\s*(\d+\.?\d*)/)||3;return{key:'prism3',geo:prism3(a,h),title:'正三棱柱 ABC-A₁B₁C₁（底边 '+a+'，高 '+h+'）'};}
    if(/ABCD.*A1|A1B1C1D1/i.test(t)){const a=num(t,/棱长[为是=]?\s*(\d+\.?\d*)/)||2;return{key:'cube',geo:cube(a),title:'正方体 ABCD-A₁B₁C₁D₁（棱长 '+a+'）',a:a};}
    const a=num(t,/棱长[为是=]?\s*(\d+\.?\d*)/)||2;return{key:'cube',geo:cube(a),title:'正方体 ABCD-A₁B₁C₁D₁（棱长 '+a+'）',a:a};
  }
  function detectQuery(t){
    let m=t.match(/二面角\s*([A-Z]\d?)\s*[-–—]\s*([A-Z]\d?[A-Z]\d?)\s*[-–—]\s*([A-Z]\d?)/);
    if(m){const e=tokPts(m[2]);if(e.length===2)return{type:'dihedral',h1:m[1],edge:e,h2:m[3]};}
    m=t.match(/(?:点)?\s*([A-Z]\d?)\s*到\s*(?:平面|底面)?\s*([A-Z]\d?(?:[A-Z]\d?)+)\s*(?:的)?\s*距离/);
    if(m){const pl=tokPts(m[2]);if(pl.length>=3)return{type:'pointplane',P:m[1],plane:pl};}
    m=t.match(/([A-Z]\d?[A-Z]\d?)\s*(?:与|和|跟)\s*(?:平面|底面)\s*([A-Z]\d?(?:[A-Z]\d?)+)\s*所?成/);
    if(m){const l=tokPts(m[1]),pl=tokPts(m[2]);if(l.length===2&&pl.length>=3)return{type:'lineplane',line:l,plane:pl};}
    m=t.match(/([A-Z]\d?[A-Z]\d?)\s*(?:与|和|跟)\s*([A-Z]\d?[A-Z]\d?)\s*所?成(?:的)?角/);
    if(!m)m=t.match(/异面直线\s*([A-Z]\d?[A-Z]\d?)\s*(?:与|和|跟)?\s*([A-Z]\d?[A-Z]\d?)/);
    if(m){const a=tokPts(m[1]),b=tokPts(m[2]);if(a.length===2&&b.length===2)return{type:'lineline',l1:a,l2:b};}
    if(/体积/.test(t))return{type:'volume'};
    return null;
  }

  function baseModel(geo){
    const span=maxSpan(geo.points);
    return {scale:1,target:threeOf(centroid(geo.points)),
      initialCamera:[span*1.9,span*1.7,span*1.9],
      points:geo.points,spheres:Object.keys(geo.points),
      edges:geo.edges.map(e=>({a:e[0],b:e[1]})),
      elements:{Axis:{type:'axes',size:span*1.25}}};
  }

  function buildLessonData(text){
    const t=normSub(text);
    const bi=detectBody(t),geo=bi.geo,P=geo.points;
    const q=detectQuery(t);
    const model=baseModel(geo);
    const span=maxSpan(P);
    const cam={x:span*1.9,y:span*1.7,z:span*1.9};
    const has=arr=>arr&&arr.every(n=>P[n]);
    let lesson,steps;

    if(q&&q.type==='lineline'&&has(q.l1)&&has(q.l2)){
      const A=P[q.l1[0]],B=P[q.l1[1]],C=P[q.l2[0]],D=P[q.l2[1]];
      const d1=V.sub(B,A),d2=V.sub(D,C),dp=V.dot(d1,d2),n1=V.len(d1),n2=V.len(d2);
      const cos=Math.abs(dp)/(n1*n2),ang=Math.acos(Math.min(1,cos))*180/Math.PI;
      const L1=q.l1.join(''),L2=q.l2.join('');
      model.elements.Line1={type:'line',a:q.l1[0],b:q.l1[1],color:'emphasis',depthTest:false};
      model.elements.Line2={type:'line',a:q.l2[0],b:q.l2[1],color:'normal',depthTest:false};
      model.elements.M1={type:'measure',a:q.l1[0],b:q.l1[1],label:L1+'='+lenLatex(n1)};
      model.elements.M2={type:'measure',a:q.l2[0],b:q.l2[1],label:L2+'='+lenLatex(n2)};
      lesson={language:'zh-CN',meta:'交互解题 · 异面直线夹角',title:bi.title+'，求异面直线 '+L1+' 与 '+L2+' 所成角',
        answerLabel:'异面直线 '+L1+' 与 '+L2+' 所成角',answerValue:'$'+deg(ang)+'$'};
      steps=[
        {title:'建立空间直角坐标系',content:'<p>以 $A$ 为原点建立坐标系，关键点：</p>$$'+[q.l1[0],q.l1[1],q.l2[0],q.l2[1]].map(n=>coordTex(n,P[n])).join(',\\ ')+'$$',highlight:['Axis'],cameraPos:cam},
        {title:'求两条直线的方向向量',content:'<p>'+vecTex(L1,d1)+'，'+vecTex(L2,d2)+'</p>',highlight:['Line1','Line2','M1','M2'],cameraPos:cam},
        {title:'用向量夹角公式',content:'<p>$$\\cos\\theta=\\frac{|\\vec{'+L1+'}\\cdot\\vec{'+L2+'}|}{|\\vec{'+L1+'}||\\vec{'+L2+'}|}=\\frac{'+numLatex(Math.abs(dp))+'}{'+lenLatex(n1)+'\\times'+lenLatex(n2)+'}='+fracTrigLatex(cos)+'$$</p>',highlight:['Line1','Line2'],cameraPos:cam},
        {title:'结论',content:'<p>异面直线 '+L1+' 与 '+L2+' 所成角 $\\theta='+deg(ang)+'$（$\\cos\\theta='+fracTrigLatex(cos)+'$）。</p>',highlight:['Line1','Line2'],cameraPos:cam}
      ];
      return {lesson:lesson,steps:steps,model:model};
    }

    if(q&&q.type==='lineplane'&&has(q.line)&&has(q.plane)){
      const A=P[q.line[0]],B=P[q.line[1]],p0=P[q.plane[0]],p1=P[q.plane[1]],p2=P[q.plane[2]];
      const nrm=V.cross(V.sub(p1,p0),V.sub(p2,p0)),d=V.sub(B,A);
      const dp=V.dot(d,nrm),nd=V.len(d),nn=V.len(nrm);
      const sin=Math.abs(dp)/(nd*nn),ang=Math.asin(Math.min(1,sin))*180/Math.PI;
      const tt=V.dot(V.sub(B,p0),nrm)/V.dot(nrm,nrm),foot=V.sub(B,V.scale(nrm,tt));
      const LN=q.line.join(''),PL=q.plane.join('');
      model.points.__F=foot;
      model.elements.Line1={type:'line',a:q.line[0],b:q.line[1],color:'emphasis',depthTest:false};
      model.elements.Plane1={type:'plane',pts:q.plane};
      model.elements.Normal1={type:'arrow',origin:q.plane[0],dir:threeOf(nrm),length:span*0.8,color:'normal'};
      model.elements.Perp={type:'line',a:q.line[1],b:'__F',color:'aux',dashed:true,depthTest:false};
      model.elements.M1={type:'measure',a:q.line[0],b:q.line[1],label:LN+'='+lenLatex(nd)};
      lesson={language:'zh-CN',meta:'交互解题 · 线面角',title:bi.title+'，求直线 '+LN+' 与平面 '+PL+' 所成角的正弦值',
        answerLabel:'直线 '+LN+' 与平面 '+PL+' 所成角的正弦值',answerValue:'$'+fracTrigLatex(sin)+'$'};
      steps=[
        {title:'建立空间直角坐标系',content:'<p>以 $A$ 为原点建系，关键点：</p>$$'+[q.line[0],q.line[1],q.plane[0],q.plane[1],q.plane[2]].map(n=>coordTex(n,P[n])).join(',\\ ')+'$$',highlight:['Axis'],cameraPos:cam},
        {title:'求平面 '+PL+' 的法向量',content:'<p>取平面内两向量叉乘得法向量 $\\vec{n}=('+nrm.map(numLatex).join(',')+')$。</p>',highlight:['Plane1','Normal1'],cameraPos:cam},
        {title:'求直线 '+LN+' 的方向向量',content:'<p>'+vecTex(LN,d)+'</p>',highlight:['Line1','M1'],cameraPos:cam},
        {title:'用向量公式求线面角',content:'<p>$$\\sin\\theta=\\frac{|\\vec{'+LN+'}\\cdot\\vec{n}|}{|\\vec{'+LN+'}||\\vec{n}|}=\\frac{'+numLatex(Math.abs(dp))+'}{'+lenLatex(nd)+'\\times'+lenLatex(nn)+'}='+fracTrigLatex(sin)+'$$</p>',highlight:['Line1','Plane1','Normal1','Perp'],cameraPos:cam},
        {title:'结论',content:'<p>直线与平面所成角 $\\theta='+deg(ang)+'$，正弦值为 $'+fracTrigLatex(sin)+'$。</p>',highlight:['Line1','Plane1','Normal1','Perp'],cameraPos:cam}
      ];
      return {lesson:lesson,steps:steps,model:model};
    }

    if(q&&q.type==='pointplane'&&P[q.P]&&has(q.plane)){
      const X=P[q.P],p0=P[q.plane[0]],p1=P[q.plane[1]],p2=P[q.plane[2]];
      const nrm=V.cross(V.sub(p1,p0),V.sub(p2,p0));
      const dist=Math.abs(V.dot(V.sub(X,p0),nrm))/V.len(nrm);
      const tt=V.dot(V.sub(X,p0),nrm)/V.dot(nrm,nrm),foot=V.sub(X,V.scale(nrm,tt));
      const PL=q.plane.join('');
      model.points.__F=foot;
      model.elements.Plane1={type:'plane',pts:q.plane};
      model.elements.Normal1={type:'arrow',origin:q.plane[0],dir:threeOf(nrm),length:span*0.8,color:'normal'};
      model.elements.Perp={type:'line',a:q.P,b:'__F',color:'emphasis',dashed:true,depthTest:false};
      model.elements.Md={type:'measure',a:q.P,b:'__F',label:'d='+lenLatex(dist)};
      lesson={language:'zh-CN',meta:'交互解题 · 点到平面距离',title:bi.title+'，求点 '+q.P+' 到平面 '+PL+' 的距离',
        answerLabel:'点 '+q.P+' 到平面 '+PL+' 的距离',answerValue:'$'+lenLatex(dist)+'$'};
      steps=[
        {title:'建立空间直角坐标系',content:'<p>关键点：</p>$$'+[q.P,q.plane[0],q.plane[1],q.plane[2]].map(n=>coordTex(n,P[n])).join(',\\ ')+'$$',highlight:['Axis'],cameraPos:cam},
        {title:'求平面 '+PL+' 的法向量',content:'<p>$\\vec{n}=('+nrm.map(numLatex).join(',')+')$</p>',highlight:['Plane1','Normal1'],cameraPos:cam},
        {title:'用点到平面距离公式',content:'<p>$$d=\\frac{|\\vec{'+q.plane[0]+q.P+'}\\cdot\\vec{n}|}{|\\vec{n}|}=\\frac{'+numLatex(Math.abs(V.dot(V.sub(X,p0),nrm)))+'}{'+lenLatex(V.len(nrm))+'}='+lenLatex(dist)+'$$</p>',highlight:['Plane1','Normal1','Perp','Md'],cameraPos:cam},
        {title:'结论',content:'<p>点 '+q.P+' 到平面 '+PL+' 的距离为 $'+lenLatex(dist)+'$。</p>',highlight:['Plane1','Perp','Md'],cameraPos:cam}
      ];
      return {lesson:lesson,steps:steps,model:model};
    }

    if(q&&q.type==='dihedral'&&P[q.h1]&&P[q.h2]&&has(q.edge)){
      const A=P[q.edge[0]],B=P[q.edge[1]],C=P[q.h1],D=P[q.h2],u=V.sub(B,A);
      const perp=Q=>{const w=V.sub(Q,A),k=V.dot(w,u)/V.dot(u,u);return V.sub(w,V.scale(u,k));};
      const v1=perp(C),v2=perp(D),cos=V.dot(v1,v2)/(V.len(v1)*V.len(v2));
      const ang=Math.acos(Math.max(-1,Math.min(1,cos)))*180/Math.PI;
      const EG=q.edge.join('');
      model.elements.Edge={type:'line',a:q.edge[0],b:q.edge[1],color:'emphasis',depthTest:false};
      model.elements.ArrU={type:'arrow',origin:q.edge[0],dir:threeOf(v1),length:span*0.7,color:'plane'};
      model.elements.ArrV={type:'arrow',origin:q.edge[0],dir:threeOf(v2),length:span*0.7,color:'normal'};
      lesson={language:'zh-CN',meta:'交互解题 · 二面角',title:bi.title+'，求二面角 '+q.h1+'-'+EG+'-'+q.h2,
        answerLabel:'二面角 '+q.h1+'-'+EG+'-'+q.h2,answerValue:'$'+deg(ang)+'$'};
      steps=[
        {title:'建立空间直角坐标系',content:'<p>关键点：</p>$$'+[q.edge[0],q.edge[1],q.h1,q.h2].map(n=>coordTex(n,P[n])).join(',\\ ')+'$$',highlight:['Axis','Edge'],cameraPos:cam},
        {title:'在两半平面内作垂直于棱的向量',content:'<p>$\\vec{v_1}=('+v1.map(numLatex).join(',')+')$，$\\vec{v_2}=('+v2.map(numLatex).join(',')+')$，均垂直于棱 '+EG+'。</p>',highlight:['Edge','ArrU','ArrV'],cameraPos:cam},
        {title:'用向量夹角求二面角',content:'<p>$$\\cos\\theta=\\frac{\\vec{v_1}\\cdot\\vec{v_2}}{|\\vec{v_1}||\\vec{v_2}|}='+numLatex(cos)+'$$</p>',highlight:['Edge','ArrU','ArrV'],cameraPos:cam},
        {title:'结论',content:'<p>二面角 '+q.h1+'-'+EG+'-'+q.h2+' $='+deg(ang)+'$。</p>',highlight:['Edge','ArrU','ArrV'],cameraPos:cam}
      ];
      return {lesson:lesson,steps:steps,model:model};
    }

    if(q&&q.type==='volume'){
      let vol=null,how='';
      if(bi.key==='cube'){const a=Math.cbrt(0)+ (P.B[0]);vol=a*a*a;how='V=a^3='+numLatex(a)+'^3';}
      else if(bi.key==='cuboid'){const lx=P.B[0],ly=P.D[1],lz=P.A1[2];vol=lx*ly*lz;how='V='+numLatex(lx)+'\\times'+numLatex(ly)+'\\times'+numLatex(lz);}
      else if(bi.key==='pyramid4'){const a=P.B[0],h=P.P[2];vol=a*a*h/3;how='V=\\frac{1}{3}S_{底}h=\\frac{1}{3}\\cdot'+numLatex(a)+'^2\\cdot'+numLatex(h);}
      else if(bi.key==='prism3'){const a=P.B[0],h=P.A1[2];vol=Math.sqrt(3)/4*a*a*h;how='V=S_{底}h=\\frac{\\sqrt{3}}{4}\\cdot'+numLatex(a)+'^2\\cdot'+numLatex(h);}
      else if(bi.key==='tetra'){vol=Math.abs(V.dot(V.sub(P.B,P.A),V.cross(V.sub(P.C,P.A),V.sub(P.D,P.A))))/6;how='V=\\frac{1}{6}|\\vec{AB}\\cdot(\\vec{AC}\\times\\vec{AD})|';}
      if(vol!=null){
        lesson={language:'zh-CN',meta:'交互解题 · 体积',title:bi.title+'，求体积',
          answerLabel:'几何体体积',answerValue:'$'+numLatex(vol)+'$'};
        steps=[
          {title:'确定几何体与尺寸',content:'<p>'+bi.title+'。</p>',highlight:['Axis'],cameraPos:cam},
          {title:'代入体积公式',content:'<p>$$'+how+'$$</p>',highlight:[],cameraPos:cam},
          {title:'结论',content:'<p>体积 $V='+numLatex(vol)+'$。</p>',highlight:[],cameraPos:cam}
        ];
        return {lesson:lesson,steps:steps,model:model};
      }
    }

    // 兜底：未识别所求，仅展示几何体
    lesson={language:'zh-CN',meta:'交互演示',title:bi.title,answerLabel:'提示',answerValue:''};
    steps=[
      {title:'三维图形',content:'<p>已生成 '+bi.title+' 的可交互三维图。左键旋转、滚轮缩放查看。</p><p>要看分步解答，请把所求写清楚，例如：<br>“正方体ABCD-A1B1C1D1棱长2，求直线A1C与底面ABCD所成角”<br>“…求异面直线AC与BD1所成角”<br>“…求点B1到平面ABCD的距离”</p>',highlight:['Axis'],cameraPos:cam}
    ];
    return {lesson:lesson,steps:steps,model:model};
  }

  /* ---- 浏览器引导：读 URL ?q= 构造数据并注入数据岛；点解答即带参重载 ---- */
  if(typeof document!=='undefined'){
    const DEF='正方体ABCD-A1B1C1D1棱长为2，求直线A1C与底面ABCD所成角';
    let qstr=DEF;
    try{const u=new URLSearchParams(location.search);if(u.get('q'))qstr=u.get('q');}catch(e){}
    window.__geoQ=qstr;
    try{
      const DATA=buildLessonData(qstr);
      const island=document.getElementById('lesson-data');
      if(island)island.textContent=JSON.stringify(DATA);
    }catch(e){console.error('buildLessonData 失败',e);}
    window.geoSolve=function(){
      const i=document.getElementById('geo-input');
      const v=i?i.value.trim():'';
      if(v)location.search='?q='+encodeURIComponent(v);
    };
    document.addEventListener('DOMContentLoaded',function(){
      const i=document.getElementById('geo-input');if(i)i.value=window.__geoQ;
      const b=document.getElementById('geo-btn');if(b)b.addEventListener('click',window.geoSolve);
      if(i)i.addEventListener('keydown',function(e){if(e.key==='Enter')window.geoSolve();});
    });
  }
  if(typeof module!=='undefined')module.exports={buildLessonData:buildLessonData,V:V,lenLatex:lenLatex,fracTrigLatex:fracTrigLatex};
