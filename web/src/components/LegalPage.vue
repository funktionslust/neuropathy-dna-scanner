<script setup lang="ts">
defineEmits<{
  close: [];
}>();
</script>

<template>
  <section class="card legal">
    <button class="close-btn" @click="$emit('close')">&times;</button>

    <h2>Legal Notice</h2>
    <p>
      <strong>Funktionslust GmbH</strong><br />
      Brudermuehlstrasse 36<br />
      81371 Munich, Germany
    </p>
    <p>
      Commercial register: HRB 241736, Munich District Court<br />
      Managing director: Wolfgang Stark<br />
      VAT ID: DE318856817
    </p>
    <p>
      Responsible for content:<br />
      Wolfgang Stark, Brudermuehlstrasse 36, 81371 Munich
    </p>
    <p>
      Questions, bug reports, or partnership enquiries:
      <a href="mailto:info@funktionslust.digital">info@funktionslust.digital</a>
    </p>

    <h2>Why this tool exists</h2>
    <p>
      Neuropathy DNA Scanner is a technical showcase by
      <a href="https://funktionslust.com" target="_blank" rel="noopener">Funktionslust GmbH</a>,
      a software development studio based in Munich. We build custom web and
      application software for clients across Europe. This tool demonstrates
      what a modern browser can do: parsing 50+ GB CRAM files, running a
      CIGAR-aware pileup genotyper, and decoding a reference-based
      compression format entirely client-side, without a backend.
    </p>
    <p>
      We built it because it needed to exist. There was no free, private,
      locally-runnable screener for inherited neuropathies. Publishing it as
      open-source lets anyone with a consumer WGS file check for the most
      common forms of inherited neuropathy without uploading their genome
      to a third party. If you have a technical challenge of your own, we
      are happy to help:
      <a href="mailto:info@funktionslust.digital">info@funktionslust.digital</a>.
    </p>

    <h2>Privacy</h2>

    <h3>What this tool does with your data</h3>
    <p>
      <strong>Nothing leaves your computer.</strong> Neuropathy DNA Scanner runs
      entirely in your browser using WebAssembly. Your genome files (BAM, CRAM,
      FASTA) are read locally by JavaScript and never transmitted to any server.
      There is no upload, no cloud processing, no backend.
    </p>

    <h3>No personal data collected</h3>
    <p>
      This tool does not collect, store, or process any personal data. There are:
    </p>
    <ul>
      <li>No user accounts or registration</li>
      <li>No cookies</li>
      <li>No analytics or tracking</li>
      <li>No server-side logging of usage</li>
      <li>No external API calls during analysis</li>
    </ul>
    <p>
      The only local storage used is a single flag remembering that you accepted
      the disclaimer, so you don't have to click it again. This is stored in
      your browser's localStorage and never sent anywhere. You can clear it
      anytime via your browser settings.
    </p>
    <p>
      The static website files (HTML, CSS, JavaScript, WebAssembly) are served
      from a web server. Standard server access logs may record your IP address
      and browser type when loading the page, as is normal for any website visit.
      Your genome data and analysis results are never sent to any server.
    </p>

    <h3>Third-party data sources</h3>
    <p>
      The variant catalog embedded in this tool was generated from publicly
      available data:
    </p>
    <ul>
      <li>
        <a href="https://www.ncbi.nlm.nih.gov/clinvar/" target="_blank" rel="noopener">ClinVar</a>
        (NCBI) - variant clinical significance
      </li>
      <li>
        <a href="https://panelapp.genomicsengland.co.uk/" target="_blank" rel="noopener">PanelApp</a>
        (Genomics England) - gene panel curation
      </li>
    </ul>
    <p>
      These databases are queried at build time only. No connections to external
      services are made when you use the tool.
    </p>

    <h3>Your rights</h3>
    <p>
      Since no personal data is collected or processed by this tool, GDPR data
      subject rights (access, rectification, erasure, portability) do not apply
      to its use. For questions about the static website hosting, contact
      Funktionslust GmbH at the address above.
    </p>

    <h3>Open source</h3>
    <p>
      The complete source code is available under the MIT licence. You can verify
      that no data is transmitted by inspecting the code or monitoring network
      traffic in your browser's developer tools.
    </p>

    <div class="action">
      <button class="btn-secondary" @click="$emit('close')">Close</button>
    </div>
  </section>
</template>

<style scoped>
.legal {
  max-width: 700px;
  margin: 0 auto;
}
.legal h2 {
  margin-top: 1.5rem;
}
.legal h2:first-of-type {
  margin-top: 0;
}
.legal h3 {
  font-size: 0.95rem;
  margin: 1rem 0 0.3rem;
}
.legal p,
.legal li {
  font-size: 0.85rem;
  line-height: 1.5;
}
.legal ul {
  margin: 0.3rem 0 0.5rem 1.2rem;
}
.close-btn {
  float: right;
  background: none;
  border: none;
  font-size: 1.5rem;
  cursor: pointer;
  color: var(--text-light);
  line-height: 1;
}
.action {
  text-align: center;
  margin-top: 1.5rem;
}
</style>
