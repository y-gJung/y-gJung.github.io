(() => {
  'use strict';
  // 24-particle Al-Salam–Carlitz ensemble on the two-sided q-lattice.
  // All particle glyphs are flat solid fills: no gradients, gloss, highlights or shadows.
  // The visible UI is deliberately label-free, but the maths remains explicit here.
  const A = -2/3;
  const N = 24;
  const $ = id => document.getElementById(id);
  const root = document.querySelector('.q-lattice');
  const visual = $('latticeVisual');
  const canvas = $('latticeCanvas');
  const ctx = canvas.getContext('2d');
  const slider = $('latticeLambda');
  const playButton = $('latticePlay');
  const resampleButton = $('latticeResample');
  const announcer = $('latticeAnnounce');
  const lMin = Number(slider.min), lMax = Number(slider.max);
  const lambda1 = Math.log(1 - A);
  const lambda2 = Math.log(1 - A) - Math.log(-A);
  const prefersReducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  const icons = {
    pause: '<path d="M9 5v14M15 5v14"/>',
    play: '<path d="m8 5 11 7-11 7V5Z"/>'
  };
  let playing = !prefersReducedMotion;
  let inView = true, model = null, state = [], pinned = null, hovered = -1;
  let from = [], target = [], animated = [];
  let animStarted = 0, lastSampled = 0, lambda = Number(slider.value);
  // The slider is the target; only the displayed density parameter is eased.
  // This leaves the exact limiting-density formula unchanged at every intermediate lambda.
  let densityLambda = lambda, lastDensityFrame = 0;
  // A track click moves the thumb at a CONSTANT rate in lambda units per second.
  // The density has its own display parameter and can ease near edges independently.
  // Direct thumb drags and keyboard/touch remain native and responsive.
  let clickGoal = null, trackPointerId = null, lastTrackFrame = 0;
  let lastTrackRebuild = 0, rebuiltAtLambda = lambda;
  const TRACK_THUMB_PX = 12;

  // Independent tuning: TRACK_SPEED affects ONLY the slider thumb after track clicks.
  // Density easing never changes the actual slider position or underlying model value.
  const TRACK_SPEED = 0.90;        // lambda/second; no slowdown near the transitions
  const DENSITY_FOLLOW_MS = 220;   // smoothing of the displayed limiting density
  const DENSITY_MAX_SPEED = 1.30;  // maximum visual lambda units/second away from edges
  const DENSITY_EDGE_SPEED = 0.46; // lower cap for the CURVE only, not the thumb
  const DENSITY_EDGE_WIDTH = 0.085;
  const DENSITY_SCALE_MS = 420;
  function maxDensitySpeed(value) {
    const distance = Math.min(Math.abs(value - lambda1), Math.abs(value - lambda2));
    const t = Math.max(0, Math.min(1, distance / DENSITY_EDGE_WIDTH));
    const smooth = t*t*(3 - 2*t);
    return DENSITY_EDGE_SPEED + (DENSITY_MAX_SPEED - DENSITY_EDGE_SPEED)*smooth;
  }
  function clamp01(x) { return Math.max(0, Math.min(1, x)); }
  function trackValueFromEvent(event) {
    const rect = slider.getBoundingClientRect();
    const fraction = clamp01((event.clientX - rect.left - TRACK_THUMB_PX/2) /
      Math.max(1, rect.width - TRACK_THUMB_PX));
    return lMin + fraction*(lMax-lMin);
  }
  function scheduleRebuild() {
    clearTimeout(pendingTimer);
    pendingTimer = setTimeout(rebuild, 145);
  }
  function advanceTrack(now) {
    const dt = lastTrackFrame ? Math.min(50, Math.max(0, now-lastTrackFrame)) : 16.7;
    lastTrackFrame = now;
    if (clickGoal === null) return;
    const distance = clickGoal - lambda;
    // Linear movement in lambda gives a uniform pixel speed along the range.
    const move = Math.sign(distance)*Math.min(Math.abs(distance), TRACK_SPEED*dt/1000);
    slider.value = String(lambda + move); // Native range quantizes to the step.
    lambda = Number(slider.value);        // Keep the model value aligned with the thumb.
    updateSlider();
    if (Math.abs(clickGoal-lambda) <= Number(slider.step)/2 + 0.000001) {
      slider.value = String(clickGoal);
      lambda = Number(slider.value);
      updateSlider();
      clickGoal = null;
      scheduleRebuild();
    } else if (now-lastTrackRebuild > 750 && Math.abs(lambda-rebuiltAtLambda) > .095) {
      // Let the particles follow on the way, rather than jumping only at the end.
      // Keep sampler work sparse to protect mobile performance.
      rebuild(N*220);lastTrackRebuild=now;
    }
  }
  let width = 800, height = 320, pendingTimer = 0, densityPath = [];
  let densityScale = 1, targetDensityScale = 1;
  let theme = {};

  function css(name) { return getComputedStyle(root).getPropertyValue('--' + name).trim() || getComputedStyle(document.documentElement).getPropertyValue('--' + name).trim(); }
  function readTheme() {
    theme = Object.fromEntries(['axis','zero-guide','site','particle','density','density-fill','accent'].map(key=>[key,css(key)]));
  }
  function updateSlider() {
    const fraction = (lambda - lMin)/(lMax - lMin);
    slider.style.setProperty('--progress', `${(fraction*100).toFixed(2)}%`);
    slider.setAttribute('aria-valuetext', `lambda ${lambda.toFixed(2)}`);
    for (const [id, critical] of [['latticeCrit1',lambda1], ['latticeCrit2',lambda2]]) {
      const mark = $(id);
      mark.style.left = `${100*(critical-lMin)/(lMax-lMin)}%`;
      mark.title = `Critical lambda ${critical.toFixed(3)}`;
    }
  }
  function setPlay(shouldPlay) {
    playing = shouldPlay;
    const label = playing ? 'Pause animation' : 'Play animation';
    playButton.setAttribute('aria-label',label);
    playButton.title = label;
    playButton.querySelector('svg').innerHTML = icons[playing ? 'pause' : 'play'];
    lastSampled = performance.now();
  }

  // Al-Salam–Carlitz weight w(x) ∝ (qx, qx/a; q)_∞.
  // The discrete Jackson q-integral contributes |x| to each lattice mass.
  function buildModel() {
    const q = Math.exp(-lambda/N);
    const K = Math.min(1600, Math.ceil(Math.log(0.0015)/Math.log(q)));
    const p = new Float64Array(K+1), aTerm = new Float64Array(K+1), invATerm = new Float64Array(K+1);
    let power = Math.pow(q,K+1), j = K+1;
    while (power>1e-13 && j<18000) {
      p[K] += Math.log1p(-power);
      aTerm[K] += Math.log1p(-A*power);
      invATerm[K] += Math.log1p(-power/A);
      power *= q; j++;
    }
    for(let k=K-1;k>=0;k--) {
      const t = Math.pow(q,k+1);
      p[k] = p[k+1] + Math.log1p(-t);
      aTerm[k] = aTerm[k+1] + Math.log1p(-A*t);
      invATerm[k] = invATerm[k+1] + Math.log1p(-t/A);
    }
    const sites=[];
    for(let k=0;k<=K;k++) {
      const factor=Math.pow(q,k);
      sites.push({ x:A*factor, k, side:-1, logWeight:Math.log(-A*factor)+p[k]+aTerm[k] });
      sites.push({ x:factor, k, side:1, logWeight:Math.log(factor)+p[k]+invATerm[k] });
    }
    const sorted=[];
    for(let k=0;k<=K;k++)sorted.push(2*k);
    for(let k=K;k>=0;k--)sorted.push(2*k+1);
    return { q, K, sites, sorted };
  }
  const interaction = (x,y) => 2*Math.log(Math.abs(x-y));
  function initialize(fixed) {
    const sites=model.sites, M=sites.length;
    const occupied=new Uint8Array(M), scores=new Float64Array(M), selected=[];
    for(let i=0;i<M;i++)scores[i]=sites[i].logWeight;
    if(fixed!==null) {
      selected.push(fixed); occupied[fixed]=1;
      for(let i=0;i<M;i++)if(i!==fixed)scores[i]+=interaction(sites[i].x,sites[fixed].x);
    }
    while(selected.length<N) {
      let max=-Infinity;
      for(let i=0;i<M;i++)if(!occupied[i]&&scores[i]>max)max=scores[i];
      let mass=0;
      for(let i=0;i<M;i++)if(!occupied[i])mass+=Math.exp(scores[i]-max);
      let r=Math.random()*mass, chosen=-1;
      for(let i=0;i<M;i++)if(!occupied[i]) {r-=Math.exp(scores[i]-max);if(r<=0){chosen=i;break;}}
      if(chosen<0){chosen=occupied.findIndex(v=>v===0);}
      occupied[chosen]=1;selected.push(chosen);
      for(let i=0;i<M;i++)if(!occupied[i])scores[i]+=interaction(sites[i].x,sites[chosen].x);
    }
    return selected;
  }
  // MCMC targets the finite truncation of the actual joint law. Samples are approximate.
  function metropolis(steps) {
    const sites=model.sites, M=sites.length, K=model.K;
    const occupied=new Set(state);
    for(let step=0;step<steps;step++) {
      const p=Math.floor(Math.random()*N), old=state[p];
      if(pinned!==null && old===pinned)continue;
      let proposed;
      if(Math.random()<.78) {
        const site=sites[old], delta=1+Math.floor(Math.pow(Math.random(),1.8)*32);
        const newK=(site.k+(Math.random()<.5?delta:-delta)+(K+1))%(K+1);
        proposed=2*newK+(site.side<0?0:1);
      }else proposed=Math.floor(Math.random()*M);
      if(proposed===old||occupied.has(proposed))continue;
      const x=sites[proposed].x, oldX=sites[old].x;
      let logRatio=sites[proposed].logWeight-sites[old].logWeight;
      for(let j=0;j<N;j++)if(j!==p) {
        const y=sites[state[j]].x;
        logRatio+=interaction(x,y)-interaction(oldX,y);
      }
      if(logRatio>=0 || Math.log(Math.random())<logRatio) {
        occupied.delete(old);occupied.add(proposed);state[p]=proposed;
      }
    }
  }
  function particles() {return state.map(i=>model.sites[i].x).sort((x,y)=>x-y);}
  function animateTo(next,previous) {
    from = previous&&previous.length===N ? previous.slice() : next.slice();
    target=next.slice();animated=from.slice();
    animStarted=performance.now();lastSampled=animStarted;
  }
  function rebuild(mcmcSteps = N*1200) {
    const previous=animated.length===N?animated.slice():null;
    model=buildModel();rebuiltAtLambda=lambda;pinned=null;hovered=-1;
    state=initialize(null);metropolis(mcmcSteps);
    makeDensity();animateTo(particles(),previous);
  }
  function resample() {
    if(!model){rebuild();return;}
    const previous=animated.length===N?animated.slice():particles();
    metropolis(N*550);
    animateTo(particles(),previous);
  }
  function pin(index) {
    pinned=pinned===index?null:index;
    const previous=animated.length===N?animated.slice():particles();
    state=initialize(pinned);metropolis(N*1300);
    animateTo(particles(),previous);
    announcer.textContent=pinned===null?'Site released':`Pinned at x = ${model.sites[pinned].x.toFixed(4)}`;
  }

  // Limiting density from Byun–Jung–Oh, Theorem 1.3 (arXiv:2507.18042), valid for -1 <= a < 0.
  function limitingDensity(x, currentLambda) {
    if(x<=A || x>=1)return 0;
    const s=Math.exp(-currentLambda), u=(1+A)*s, v=2*Math.sqrt(-A*(1-s)*s);
    const left=u-v, right=u+v;
    if(x<left || x>right) {
      if(currentLambda>lambda1 && x>right && x<1)return 1/(currentLambda*Math.abs(x));
      if(currentLambda>lambda2 && x>A && x<left)return 1/(currentLambda*Math.abs(x));
      return 0;
    }
    const xx=Math.abs(x)<1e-5?(x<0?-1e-5:1e-5):x;
    const x0=(A*A+1-xx*(A+1))/((A-1)*(A-1));
    const x1=Math.sqrt(Math.max(0,4*A*(xx-A)*(xx-1)))/((A-1)*(A-1));
    const ratio=((1-x0-x1)/(1-x0+x1))*((1-s-x0+x1)/(x0+x1-1+s));
    if(!Number.isFinite(ratio))return 0;
    return 2/(Math.PI*currentLambda*Math.abs(xx))*Math.atan(Math.sqrt(Math.max(0,ratio)));
  }
  // The displayed lambda is updated each animation frame. No debouncing for the curve.
  function makeDensity(currentLambda = densityLambda) {
    densityPath=[]; let max=0;
    for(let i=0;i<=960;i++) {
      const x=A+(0.00001+(1-0.00002)*i/960)*(1-A);
      const d=Math.max(0,limitingDensity(x,currentLambda));
      if(Number.isFinite(d)){
        densityPath.push([x,Math.min(80,d)]);
        max=Math.max(max,Math.min(80,d));
      }
    }
    targetDensityScale=(height*.48)/Math.max(.1,max);
  }
  function advanceDensity(now) {
    const dt = lastDensityFrame ? Math.max(0, Math.min(50, now-lastDensityFrame)) : 16.7;
    lastDensityFrame = now;
    if (prefersReducedMotion) {
      if (densityLambda !== lambda) {
        densityLambda = lambda;
        makeDensity();
      }
      densityScale = targetDensityScale;
      return;
    }
    // Curve-only easing: the thumb remains independent and moves uniformly.
    // Each intermediate curve is still evaluated from the limiting-density formula.
    const follow = 1-Math.exp(-dt/DENSITY_FOLLOW_MS);
    if (Math.abs(densityLambda-lambda) > .00002) {
      const desired = (lambda-densityLambda)*follow;
      const limit = maxDensitySpeed(densityLambda)*dt/1000;
      densityLambda += Math.sign(desired)*Math.min(Math.abs(desired), limit);
      if (Math.abs(densityLambda-lambda) <= .00002) densityLambda = lambda;
      makeDensity();
    }
    // Prevent the canvas auto-scale from jolting near a hard/soft edge transition.
    densityScale += (targetDensityScale-densityScale)*(1-Math.exp(-dt/DENSITY_SCALE_MS));
  }
  const mapX=x=>30+(x-A)/(1-A)*(width-60);
  function resize() {
    const rect=canvas.getBoundingClientRect();
    width=Math.max(260,rect.width);height=Math.max(190,rect.height);
    const dpr=Math.min(2,window.devicePixelRatio||1);
    canvas.width=Math.round(width*dpr);
    canvas.height=Math.round(height*dpr);
    ctx.setTransform(dpr,0,0,dpr,0,0);
    makeDensity();
    densityScale = targetDensityScale; // Resize should not animate from the old canvas height.
  }
  function dot(x,y,r,color) {ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ctx.fillStyle=color;ctx.fill();}
  function draw(now) {
    if(inView&&!document.hidden){
      advanceTrack(now);
      advanceDensity(now);
      const y0=Math.round(height*.72)+.5; // Every lattice site and particle shares exactly this y-coordinate.
      ctx.clearRect(0,0,width,height);
      ctx.beginPath();ctx.moveTo(30,y0);ctx.lineTo(width-30,y0);
      ctx.strokeStyle=theme.axis;ctx.lineWidth=1;ctx.stroke();
      const zero=mapX(0);
      ctx.setLineDash([2,6]);ctx.beginPath();ctx.moveTo(zero,31);ctx.lineTo(zero,height-28);
      ctx.strokeStyle=theme['zero-guide'];ctx.stroke();ctx.setLineDash([]);
      if(model) {
        // Always draw the limiting density and visible q-lattice sites.
        ctx.beginPath();ctx.moveTo(mapX(A),y0);
        for(const [x,d] of densityPath)ctx.lineTo(mapX(x),y0-d*densityScale);
        ctx.lineTo(mapX(1),y0);ctx.closePath();ctx.fillStyle=theme['density-fill'];ctx.fill();
        ctx.beginPath();
        densityPath.forEach(([x,d],i)=>i?ctx.lineTo(mapX(x),y0-d*densityScale):ctx.moveTo(mapX(x),y0-d*densityScale));
        ctx.strokeStyle=theme.density;ctx.lineWidth=1.5;ctx.stroke();
        let previousX=-Infinity;
        for(const index of model.sorted){
          const x=mapX(model.sites[index].x);
          // Drawing every mathematically distinct site would merge into a solid blob near zero.
          if(x-previousX<2.5&&index!==pinned&&index!==hovered)continue;
          previousX=x;dot(x,y0,1.5,theme.site);
        }
        if(hovered>=0) {
          const x=mapX(model.sites[hovered].x);
          ctx.beginPath();ctx.arc(x,y0,7,0,Math.PI*2);ctx.strokeStyle=theme.accent;ctx.lineWidth=1;ctx.stroke();
        }
        if(pinned!==null) {
          ctx.beginPath();ctx.arc(mapX(model.sites[pinned].x),y0,10.5,0,2*Math.PI);
          ctx.strokeStyle=theme.accent;ctx.lineWidth=1;ctx.stroke();
        }
        const duration=prefersReducedMotion?1:1050;
        const t=Math.max(0,Math.min(1,(now-animStarted)/duration));
        const easing=1-Math.pow(1-t,3);
        animated=target.map((x,i)=>(from[i]??x)*(1-easing)+x*easing);
        for(const x of animated){
          const px=mapX(x);
          dot(px,y0,4.6,theme.particle);
        }
        if(playing && now-lastSampled>4800)resample();
      }
    }
    requestAnimationFrame(draw);
  }
  function nearestSite(event) {
    if(!model)return -1;
    const rect=canvas.getBoundingClientRect();
    const x=event.clientX-rect.left, y=event.clientY-rect.top;
    if(Math.abs(y-(Math.round(height*.72)+.5))>24)return -1;
    let best=10,chosen=-1;
    model.sites.forEach((site,i)=>{
      const distance=Math.abs(mapX(site.x)-x);
      if(distance<best){best=distance;chosen=i;}
    });
    return chosen;
  }
  playButton.addEventListener('click',()=>setPlay(!playing));
  resampleButton.addEventListener('click',resample);
  // Clicking the empty track normally teleports the native thumb. On desktop,
  // intercept only these jumps and animate the thumb toward the destination.
  // A grab ON the thumb, touch gestures, keyboard inputs and screen readers retain
  // standard range-input behavior.
  slider.addEventListener('pointerdown', event => {
    if (prefersReducedMotion || (event.pointerType !== 'mouse' && event.pointerType !== 'pen')) return;
    const rect = slider.getBoundingClientRect();
    const fraction = (lambda-lMin)/(lMax-lMin);
    const thumbX = rect.left + TRACK_THUMB_PX/2 + fraction*(rect.width-TRACK_THUMB_PX);
    if (Math.abs(event.clientX-thumbX) <= 12) {
      clickGoal = null; // Direct thumb drag: immediate input remains responsive.
      return;
    }
    event.preventDefault();
    slider.focus({preventScroll:true});
    trackPointerId = event.pointerId;
    slider.setPointerCapture(event.pointerId);
    const step = Number(slider.step);
    clickGoal = lMin + Math.round((trackValueFromEvent(event)-lMin)/step)*step;
    clickGoal = Math.max(lMin, Math.min(lMax, clickGoal));
    lastTrackFrame = performance.now();
    lastTrackRebuild = performance.now();
  });
  slider.addEventListener('pointermove', event => {
    if (event.pointerId === trackPointerId) {
      const step = Number(slider.step);
      clickGoal = Math.max(lMin, Math.min(lMax,
        lMin + Math.round((trackValueFromEvent(event)-lMin)/step)*step));
    }
  });
  for (const name of ['pointerup','pointercancel','lostpointercapture']) {
    slider.addEventListener(name, event => {
      if (event.pointerId === trackPointerId) trackPointerId = null;
    });
  }
  slider.addEventListener('input',()=>{
    clickGoal=null;
    lambda=Number(slider.value);updateSlider();
    // Expensive MCMC sampling is delayed until the user pauses dragging.
    scheduleRebuild();
  });
  canvas.addEventListener('pointermove',event=>{
    hovered=nearestSite(event);
    canvas.style.cursor=hovered<0?'default':'pointer';
  });
  canvas.addEventListener('pointerleave',()=>hovered=-1);
  canvas.addEventListener('click',event=>{
    const site=nearestSite(event);
    if(site>=0)pin(site); else resample();
  });
  canvas.addEventListener('keydown',event=>{
    if(event.key==='Enter'||event.key===' '){event.preventDefault();resample();}
  });
  document.addEventListener('visibilitychange',()=>{lastSampled=performance.now();});
  if('IntersectionObserver' in window) {
    new IntersectionObserver(entries=>{inView=entries[0].isIntersecting;lastSampled=performance.now();},{threshold:.01}).observe(visual);
  }
  if('ResizeObserver' in window)new ResizeObserver(resize).observe(visual);
  else window.addEventListener('resize',resize);
  readTheme();updateSlider();setPlay(playing);resize();rebuild();requestAnimationFrame(draw);
})();
