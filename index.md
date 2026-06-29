---
layout: default
title: Home
---

<section class="home-page">
  <header class="home-hero">
    <div class="home-photo-frame">
      <img src="/myphoto.jpg" alt="Portrait of Yeong-Gwang Jung">
    </div>
    <div class="home-hero-text">
      <p class="eyebrow">Welcome</p>
      <h1>Yeong-Gwang Jung</h1>
      <p class="lead">
        I am a Ph.D. student in the Department of Mathematical Sciences at Seoul National University,
        advised by <a href="https://sites.google.com/view/sungsoobyun">Sung-Soo Byun</a>.
        I am also a student member of the
        <a href="https://sites.google.com/view/snuprob/home?authuser=0">SNU Probability Group</a>.
      </p>
    </div>
  </header>

  <section class="home-section">
    <h2>Research Interests</h2>
    <div class="interest-grid">
      <div class="interest-card">
        <strong>Random Matrix Theory</strong>
        <span>Spectral statistics, limiting laws, and related probabilistic structures.</span>
      </div>
      <div class="interest-card">
        <strong>Orthogonal Polynomials</strong>
        <span>Including their <span class="math-term">\(q\)&#8209;analogues</span> and connections to special functions.</span>
      </div>
      <div class="interest-card">
        <strong>Quantum Information Theory</strong>
        <span>Mathematical aspects of quantum systems and information-theoretic problems.</span>
      </div>
      <div class="interest-card">
        <strong>Free Probability</strong>
        <span>Noncommutative probability and its interaction with random matrices.</span>
      </div>
    </div>
  </section>

  <p class="updated-line">Last updated: <span id="last-updated"></span></p>
</section>

<script>
  const updated = new Date(document.lastModified);
  document.getElementById("last-updated").textContent =
    updated.toLocaleDateString("en-US", {
      year: "numeric",
      month: "short",
      day: "numeric"
    });
</script>
