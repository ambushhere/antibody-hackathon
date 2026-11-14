#!/usr/bin/env node

/**
 * Build the three editable selection and architecture slides and PNG previews.
 *
 * Rebuild from the repository root with the supplied presentation runtime:
 *   RUNTIME_NODE=/path/to/node \
 *   RUNTIME_NODE_MODULES=/path/to/node_modules \
 *   RUNTIME_PYTHON=/path/to/python3 \
 *   PRESENTATIONS_SKILL_DIR=/path/to/presentations/skills/presentations \
 *   OUTPUT_STEM=her2-selection-and-evidence-en \
 *   "$RUNTIME_NODE" scripts/presentation/selection-slide.mjs
 *
 * The runtime module directory is linked only inside the ignored results/tmp
 * build folder, which also contains the private finalizer receipt.
 */

import fs from "node:fs/promises";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath, pathToFileURL } from "node:url";

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const workspaceDir = path.resolve(scriptDir, "../..");
const buildDir = process.env.BUILD_DIR ?? path.join(workspaceDir, "results/tmp/her2-slide-build");
const outputDir = path.join(workspaceDir, "output/selection-slide");
const outputStem = process.env.OUTPUT_STEM ?? "her2-selection-and-evidence-en";
const finalPath = path.join(outputDir, `${outputStem}.pptx`);
const candidatePath = path.join(outputDir, ".candidate-her2-pipeline.pptx");
const receiptPath = path.join(buildDir, `${path.basename(finalPath)}.validation.json`);
const skillDir = process.env.PRESENTATIONS_SKILL_DIR;
const runtimeModules = process.env.RUNTIME_NODE_MODULES ?? process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES;
const pythonExecutable = process.env.RUNTIME_PYTHON;

if (!runtimeModules || !pythonExecutable || !skillDir) {
  throw new Error("Set RUNTIME_NODE_MODULES, RUNTIME_PYTHON and PRESENTATIONS_SKILL_DIR to the supplied runtime paths.");
}

await fs.mkdir(buildDir, { recursive: true });
await fs.mkdir(outputDir, { recursive: true });
await fs.rm(candidatePath, { force: true });

const buildModules = path.join(buildDir, "node_modules");
try {
  await fs.lstat(buildModules);
} catch {
  await fs.symlink(runtimeModules, buildModules, "dir");
}

const buildRequire = createRequire(path.join(buildDir, "loader.cjs"));
const artifactToolEntry = buildRequire.resolve("@oai/artifact-tool");
const { Presentation, PresentationFile, FileBlob } = await import(pathToFileURL(artifactToolEntry).href);
const { finalizePresentation } = await import(pathToFileURL(
  path.join(skillDir, "container_tools/artifact_tool_utils.mjs"),
).href);

const W = 1280;
const H = 720;
const FONT = "Arial";
const C = {
  ink: "#10243E",
  body: "#24364B",
  muted: "#5C6D80",
  rule: "#D9E1EA",
  blue: "#2867B2",
  blueSoft: "#EDF4FC",
  green: "#267557",
  greenSoft: "#EEF6F1",
  amber: "#A45B08",
  amberSoft: "#FFF4E3",
  white: "#FFFFFF",
  transparent: "none",
};

const ppt = process.env.SOURCE_PPTX
  ? await PresentationFile.importPptx(await FileBlob.load(process.env.SOURCE_PPTX))
  : Presentation.create({ slideSize: { width: W, height: H } });

function addShape(slide, geometry, x, y, w, h, fill = C.transparent, lineFill = C.transparent, lineWidth = 0, name = undefined) {
  return slide.shapes.add({
    geometry,
    name,
    position: { left: x, top: y, width: w, height: h },
    fill,
    line: { style: "solid", fill: lineFill, width: lineWidth },
  });
}

function addText(slide, text, x, y, w, h, options = {}) {
  const shape = addShape(slide, "textbox", x, y, w, h, options.fill ?? C.transparent, options.lineFill ?? C.transparent, options.lineWidth ?? 0, options.name);
  shape.text = text;
  shape.text.style = {
    typeface: FONT,
    fontSize: options.size ?? 24,
    color: options.color ?? C.body,
    bold: options.bold ?? false,
    alignment: options.align ?? "left",
    verticalAlignment: options.valign ?? "middle",
    wrap: "square",
    autoFit: options.autoFit ?? "shrinkText",
    insets: options.insets ?? { left: 0, right: 0, top: 0, bottom: 0 },
  };
  return shape;
}

function addRule(slide, x1, y, x2, color = C.rule, thickness = 1) {
  addShape(slide, "line", x1, y, x2 - x1, 0, C.transparent, color, thickness);
}

function addArrow(slide, x, y, color) {
  addShape(slide, "rightArrow", x, y, 22, 17, color, C.transparent, 0);
}

function addStep(slide, x, y, width, label, body, accent, bodySize = 23) {
  addText(slide, label, x, y, width, 24, { size: 17, color: accent, bold: true, valign: "top" });
  addText(slide, body, x, y + 31, width, 103, { size: bodySize, color: C.ink, valign: "top" });
}

function addHeader(slide, title, subtitle) {
  slide.background.fill = C.white;
  addText(slide, title, 64, 37, 1152, 57, { size: 45, color: C.ink, bold: true, valign: "top" });
  addText(slide, subtitle, 66, 105, 1148, 37, { size: 23, color: C.muted, valign: "top" });
  addRule(slide, 64, 157, 1216, C.rule, 1.4);
}

if (!process.env.SOURCE_PPTX) {
const slide = ppt.slides.add();
addHeader(slide, "What we selected and why", "Two routes to the same HER2 binding site, with different evidence for selection.");
addText(slide, "10 optimization variants", 64, 190, 465, 45, {size:32,bold:true,color:C.blue});
addText(slide, "Refine an existing antibody", 64, 246, 465, 35, {size:25,bold:true});
addText(slide, "Keep the trastuzumab scaffold.\nChange its CDRH3 binding loop.",64,291,465,72,{size:24,valign:"top"});
addText(slide, "Choose using FACS / SPR predictions\nand sequence diversity.",64,377,465,72,{size:24,valign:"top"});
addText(slide, "Why: measured data from related variants\nhelp guide the choice.",64,462,465,70,{size:23,color:C.blue,valign:"top"});
addText(slide, "20 de novo variants", 655,190,560,45,{size:32,bold:true,color:C.green});
addText(slide, "Explore new antibody designs",655,246,560,35,{size:25,bold:true});
addText(slide, "Generate new heavy/light chain pairs\nwith RFD4-Proteina.",655,291,560,72,{size:24,valign:"top"});
addText(slide, "Choose fewer motif warnings, greater\nsequence diversity and novelty.",655,377,560,72,{size:24,valign:"top"});
addText(slide, "Why: explore alternative designs.\nAffinity ranking remains unvalidated.",655,462,560,70,{size:23,color:C.green,valign:"top"});
addShape(slide,"downArrow",290,535,24,28,C.blue);
addShape(slide,"downArrow",920,535,24,28,C.green);
addText(slide,"30 ranked candidates • All pass installed abguard",64,574,1152,38,{size:29,bold:true,align:"center"});
addText(slide,"First 4 slots: 3 optimization + 1 de novo, prioritizing the route with measured-data support.",64,618,1152,30,{size:22,align:"center"});
addText(slide,"Laboratory tests must establish binding, epitope preservation and thermal stability.",64,670,1152,26,{size:20,color:C.muted,align:"center"});
slide.speakerNotes.text = "We refine a known antibody and explore new designs in parallel. Optimization has measured data from related trastuzumab variants. De novo selection emphasizes eligibility and diversity because we have no validated affinity ranker for redesigned frameworks. Abguard checks sequence eligibility. All biological gains remain unmeasured. FACS is fluorescence-activated cell sorting; SPR is surface plasmon resonance. Sources: docs/EXECUTION_EN.md; docs/research/ranked-draft-policy.md; results/abguard-ranked/receipt.json.";

const dataSlide = ppt.slides.add();
addHeader(dataSlide, "Tools and data behind our choices", "Measured data guide optimization. Generation and eligibility checks support both tracks.");
addText(dataSlide,"Measured datasets",64,188,530,42,{size:31,bold:true,color:C.blue});
addText(dataSlide,"524,346 FACS variants",64,251,530,34,{size:26,bold:true});
addText(dataSlide,"Sorting-signal labels for CDRH3 classification",64,292,550,42,{size:22});
addText(dataSlide,"158 variants in the SPR evaluation",64,354,550,34,{size:26,bold:true});
addText(dataSlide,"From 190 quantified rows in a 768-row source.\nGrouped validation tests affinity ranking.",64,395,550,64,{size:22,valign:"top"});
addText(dataSlide,"36 external SPR controls",64,482,530,34,{size:26,bold:true});
addText(dataSlide,"25 numerical KD values for a limited transfer check",64,523,550,58,{size:22,valign:"top"});
addText(dataSlide,"Tools and their roles",655,188,560,42,{size:31,bold:true,color:C.green});
addText(dataSlide,"FACS classifiers + SPR ridge models",655,251,560,34,{size:26,bold:true});
addText(dataSlide,"Score optimization proposals on separate signals",655,292,560,42,{size:22});
addText(dataSlide,"Apheris RFD4-Proteina",655,354,560,34,{size:26,bold:true});
addText(dataSlide,"Generate 64 new heavy/light chain pairs",655,395,560,42,{size:22});
addText(dataSlide,"Abguard + immunum",655,482,560,34,{size:26,bold:true});
addText(dataSlide,"Check eligibility and antibody numbering.\nBundled novelty reference: 610,109 CDRH3s.",655,523,560,64,{size:22,valign:"top"});
addRule(dataSlide,64,610,1216);
addText(dataSlide,"Exploratory: Transformers and 100,000 synthetic labels did not replace the SPR reference.",64,628,1152,29,{size:21});
addText(dataSlide,"Structural controls did not justify candidate ranking. Thermal improvement remains unmeasured.",64,668,1152,28,{size:21,color:C.muted});
dataSlide.speakerNotes.text="FACS class is a sorting signal, not numeric affinity. The 158 SPR variants form a filtered cohort, with a nested distant subset of 75. External SPR controls are a separate limited transfer check. Synthetic labels add no measured evidence. RFD4 generated 64 pairs; the final ranked draft contains 20 selected de novo candidates. Abguard is eligibility only. Additional investigations included AbLang2, OpenDDE, Protenix, RF3 and ThermoMPNN-D; they did not justify candidate promotion. Sources: docs/EXPERIMENT_INDEX.md; docs/EXECUTION_EN.md; docs/research/neural-results.md; docs/research/abguard-validation.md.";

}
if (ppt.slides.items.length !== 2) throw new Error("Architecture append expects a two-slide source");
const arch = ppt.slides.add();
addHeader(arch, "The architecture used for optimization", "Two separate sequence scores guide selection. The SPR scorer is a small additive model.");
addText(arch,"CDRH3 window",64,198,260,35,{size:27,bold:true,color:C.blue});
addText(arch,"10 amino acids\nW G G D G F Y A M D",64,252,285,70,{size:24,valign:"top"});
addArrow(arch,356,263,C.blue);
addText(arch,"Position encoding",410,198,300,35,{size:27,bold:true,color:C.blue});
addText(arch,"10 positions × 20 letters\n= 200 one-hot features",410,252,310,70,{size:24,valign:"top"});
addArrow(arch,755,263,C.blue);
addText(arch,"SPR ridge regression",817,198,400,35,{size:27,bold:true,color:C.blue});
addText(arch,"Weighted sum + intercept\nPredicts −log10(KD in nM)",817,252,400,70,{size:24,valign:"top"});
addText(arch,"Training: 158 measured variants. Ridge shrinks weights to limit overfitting.",64,354,1152,36,{size:25,bold:true});
addText(arch,"Higher score predicts tighter binding. Candidate KD remains uncalibrated.",64,402,1152,34,{size:23,color:C.muted});
addRule(arch,64,460,1216);
addText(arch,"Parallel FACS score",64,485,300,34,{size:26,bold:true,color:C.green});
addText(arch,"Position + amino-acid pair features\nLogistic regression estimates the high class",64,530,575,70,{size:23,valign:"top"});
addArrow(arch,650,551,C.green);
addText(arch,"Selection",710,485,500,34,{size:26,bold:true,color:C.green});
addText(arch,"Balance FACS and SPR scores\nRetain sequence diversity and eligibility",710,530,505,70,{size:23,valign:"top"});
addText(arch,"Transformer: trained and benchmarked, but excluded from final selection. This scorer does not rank de novo designs.",64,656,1152,46,{size:20,color:C.muted});
arch.speakerNotes.text="The final optimization shortlist used position-one-hot ridge, with 200 coefficients and an intercept, alpha 10 selected by grouped tuning. Each sequence activates ten positional coefficients. The formula adds position effects independently; pairwise ridge was benchmarked later and did not replace the frozen shortlist. FACS uses a separate pairwise logistic model. The WT sequence is illustrative, not a submitted novel candidate. The large experimental Transformer used 10 residue tokens plus CLS, learned 256-dimensional embeddings and positions, six encoder blocks, eight attention heads and a 1024-wide feed-forward block. It predicted three FACS classes and was not promoted to final selection. Sources: scripts/rank_spr_candidates.py; src/her2_benchmark.py; results/spr-expanded-ranking/report.md; scripts/train_her2_neural.py; docs/research/neural-results.md.";

if (ppt.slides.items.length !== 3) {
  throw new Error(`Expected exactly 3 slides, got ${ppt.slides.items.length}`);
}

const expectedSlideSizeEmu = "12192000,6858000";
await (await PresentationFile.exportPptx(ppt)).save(candidatePath);

const requirements = {
  explicitTotalSlideCount: 3,
  requiredNativeTableOwnerSlides: [],
  requiredNativeChartOwnerSlides: [],
  fontPolicy: { basis: "design", families: [FONT] },
};

try {
  await finalizePresentation({
    ...requirements,
    workspaceDir,
    candidatePath,
    finalPath,
    pythonExecutable,
    integrityValidatorPath: path.join(skillDir, "container_tools/inspect_presentation_package_integrity.py"),
    layoutValidatorPath: path.join(skillDir, "container_tools/inspect_presentation_layout_geometry.py"),
    layoutArgs: [
      "--expected-slide-size-emu", expectedSlideSizeEmu,
      "--validate-heading-fit",
    ],
    fontPolicy: requirements.fontPolicy,
    verifyArtifactToolImport: true,
    receiptPath,
  });
} finally {
  await fs.rm(candidatePath, { force: true });
  await fs.rm(`${candidatePath}.inspect.ndjson`, { force: true });
}

// Render from the finalized PPTX, so these previews represent the deliverable.
const finalDeck = await PresentationFile.importPptx(await FileBlob.load(finalPath));
if (finalDeck.slides.items.length !== 3) {
  throw new Error(`Final PPTX has ${finalDeck.slides.items.length} slides, expected 3.`);
}
for (const [index, slide] of finalDeck.slides.items.entries()) {
  const image = await finalDeck.export({ slide, format: "png", scale: 1 });
  await fs.writeFile(path.join(outputDir, `slide-${String(index + 1).padStart(2, "0")}.png`), new Uint8Array(await image.arrayBuffer()));
}

console.log(JSON.stringify({
  slideCount: finalDeck.slides.items.length,
  finalPath,
  previews: [1, 2, 3].map(n => path.join(outputDir, `slide-0${n}.png`)),
  validationReceipt: receiptPath,
}, null, 2));
