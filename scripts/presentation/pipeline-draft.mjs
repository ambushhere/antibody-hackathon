#!/usr/bin/env node

/**
 * Build the two editable pipeline slides and PNG previews.
 *
 * Rebuild from the repository root with the supplied presentation runtime:
 *   RUNTIME_NODE=/path/to/node \
 *   RUNTIME_NODE_MODULES=/path/to/node_modules \
 *   RUNTIME_PYTHON=/path/to/python3 \
 *   PRESENTATIONS_SKILL_DIR=/path/to/presentations/skills/presentations \
 *   OUTPUT_STEM=her2-hackathon-pipeline-draft-v6 \
 *   "$RUNTIME_NODE" scripts/presentation/pipeline-draft.mjs
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
const outputDir = path.join(workspaceDir, "output/pipeline-draft");
const outputStem = process.env.OUTPUT_STEM ?? "her2-hackathon-pipeline-draft-v6";
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

const ppt = Presentation.create({ slideSize: { width: W, height: H } });

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

// Slide 1: both generation tracks converge on one research pool.
{
  const slide = ppt.slides.add();
  addHeader(
    slide,
    "Generating HER2 binders",
    "Shared goal: retain trastuzumab’s HER2 site and improve thermal stability, to be measured later.",
  );

  const laneX = 64;
  const laneW = 1152;
  addShape(slide, "rect", laneX, 178, laneW, 170, C.blueSoft);
  addShape(slide, "rect", laneX, 365, laneW, 170, C.greenSoft);
  addShape(slide, "rect", laneX, 178, 5, 170, C.blue);
  addShape(slide, "rect", laneX, 365, 5, 170, C.green);

  addText(slide, "OPTIMIZATION", 85, 201, 145, 30, { size: 18, bold: true, color: C.blue, valign: "top" });
  addText(slide, "Trastuzumab\nscaffold", 85, 239, 145, 67, { size: 23, bold: true, color: C.ink, valign: "top" });
  addStep(slide, 254, 197, 182, "CHANGE 10-AA CDRH3", "Statistical proposals; framework fixed", C.blue, 23);
  addArrow(slide, 444, 257, C.blue);
  addStep(slide, 479, 197, 213, "RANK", "FACS signal plus direct-SPR ridge score", C.blue, 23);
  addArrow(slide, 700, 257, C.blue);
  addStep(slide, 735, 197, 210, "SCREEN", "Local sequence filters select 10 provisional variants", C.blue, 22);
  addArrow(slide, 956, 257, C.blue);
  addText(slide, "10", 993, 217, 91, 64, { size: 49, bold: true, color: C.blue, align: "center" });
  addText(slide, "optimization\npicks", 1080, 223, 120, 55, { size: 21, bold: true, color: C.ink, valign: "middle" });

  addText(slide, "DE NOVO", 85, 388, 145, 30, { size: 18, bold: true, color: C.green, valign: "top" });
  addText(slide, "Paired VH/VL\ndesign", 85, 426, 145, 67, { size: 23, bold: true, color: C.ink, valign: "top" });
  addStep(slide, 254, 384, 182, "GENERATE", "RFD4-Proteina: antibody and HER2-novel checkpoints", C.green, 21);
  addArrow(slide, 444, 444, C.green);
  addStep(slide, 479, 384, 213, "REDESIGN", "Paired VH/VL toward the trastuzumab-site context", C.green, 22);
  addArrow(slide, 700, 444, C.green);
  addStep(slide, 735, 384, 210, "AUDIT", "64 generated → 43 unique local sequence passes", C.green, 22);
  addArrow(slide, 956, 444, C.green);
  addText(slide, "43", 993, 404, 91, 64, { size: 49, bold: true, color: C.green, align: "center" });
  addText(slide, "de novo\npairs", 1080, 410, 120, 55, { size: 21, bold: true, color: C.ink, valign: "middle" });

  addRule(slide, 64, 555, 1216, C.rule, 1.2);
  addText(slide, "52", 68, 574, 83, 59, { size: 46, bold: true, color: C.ink });
  addText(slide, "pass installed abguard\n11 novel + 41 de_novo", 158, 578, 426, 56, { size: 22, color: C.ink, bold: true });
  addText(slide, "53 screened against bundled reference. Final list ≤20; laboratory outcomes unmeasured.", 617, 584, 585, 48, { size: 18, color: C.amber, valign: "middle" });

  slide.speakerNotes.text = [
    "Питч: в оптимизации сохраняем каркас VH/VL трастузумаба и меняем десятиостаточное окно CDRH3. Сигнал FACS и ridge, обученный на измерениях SPR, дают предварительный порядок для десяти вариантов.",
    "De novo ветка использует RFD4-Proteina с antibody и HER2-novel checkpoint. Из 64 созданных пар после локальной проверки осталось 43 уникальные последовательности. Это происхождение дизайнов, не присвоенный трек. Abguard 0.1.0 с базой 610109 CDRH3 допускает 52/53 (11 novel, 41 de_novo), один дизайн отклонён. Допуск не доказывает связывание.",
    "Sources: docs/research/abguard-validation.md; docs/research/denovo-scale-results.md; results/denovo-scale-audit/summary.json; results/spr-expanded-ranking/report.md; configs/her2-candidate-screen/manifest.json; docs/research/candidate-screen-results.md.",
  ].join("\n\n");
}

// Slide 2: separate pre-lab evidence from the direct laboratory outcomes.
{
  const slide = ppt.slides.add();
  addHeader(
    slide,
    "Evaluation before the laboratory",
    "Three separate evidence axes; direct assays remain the final test.",
  );

  const cols = [64, 454, 844];
  const colW = 346;
  addRule(slide, 432, 184, 432, C.rule, 1.2);
  addRule(slide, 822, 184, 822, C.rule, 1.2);

  addText(slide, "01  SEQUENCE + SPR", cols[0], 181, colW, 28, { size: 18, bold: true, color: C.blue, valign: "top" });
  addText(slide, "Does sequence ranking transfer?", cols[0], 222, colW, 58, { size: 26, bold: true, color: C.ink, valign: "top" });
  addText(slide, "FACS: frozen no-overlap split\nGrouped SPR cross-validation\n158 measured KD variants", cols[0], 291, colW, 92, { size: 23, color: C.body, valign: "top" });
  addText(slide, "External SPR: 36 eligible controls; advantage uncertain.", cols[0], 390, colW, 48, { size: 18, color: C.muted, valign: "top" });

  addText(slide, "02  STRUCTURE + EPITOPE", cols[1], 181, colW, 28, { size: 18, bold: true, color: C.green, valign: "top" });
  addText(slide, "Controls before candidates", cols[1], 222, colW, 58, { size: 26, bold: true, color: C.ink, valign: "top" });
  addText(slide, "Protenix v2 · HER2-only template\nWT + strong / weak SPR controls\nThree matched seeds", cols[1], 291, colW, 92, { size: 23, color: C.body, valign: "top" });
  addText(slide, "Pair ipTM, contacts, PAE, RMSD and clashes; cross-model check planned.", cols[1], 390, colW, 48, { size: 18, color: C.muted, valign: "top" });

  addText(slide, "03  THERMAL STABILITY", cols[2], 181, colW, 28, { size: 18, bold: true, color: C.amber, valign: "top" });
  addText(slide, "Keep Tm unresolved", cols[2], 222, colW, 58, { size: 26, bold: true, color: C.ink, valign: "top" });
  addText(slide, "Candidate thermal evidence is missing.\nNo validated Tm score enters ranking.\nMeasure construct-matched Tm directly.\nDo not infer Tm from ipTM.", cols[2], 291, colW, 132, { size: 23, color: C.body, valign: "top" });

  addRule(slide, 64, 452, 1216, C.rule, 1.2);
  addText(slide, "STRUCTURE GATE", 66, 468, 167, 28, { size: 18, bold: true, color: C.amber, valign: "top" });
  addText(slide, "WT gate: no-MSA 1/3; target-MSA 0/3; need 2/3 → HOLD", 244, 462, 946, 38, { size: 25, bold: true, color: C.amber, valign: "top" });
  addText(slide, "Both candidate co-folding screens were withheld; no candidate poses were audited.", 66, 507, 1146, 31, { size: 20, color: C.body, valign: "top" });
  addText(slide, "Per-seed pass: ≥13/19 contacts, target RMSD ≤2 Å, interface RMSD ≤4 Å and zero clashes.", 66, 540, 1146, 27, { size: 17, color: C.muted, valign: "top" });

  addRule(slide, 64, 586, 1216, C.rule, 1.2);
  addText(slide, "LAB MEASURES", 66, 603, 166, 27, { size: 18, bold: true, color: C.ink, valign: "top" });
  addText(slide, "SPR KD", 253, 596, 204, 45, { size: 25, bold: true, color: C.blue, align: "center" });
  addText(slide, "Trastuzumab competition", 481, 596, 369, 45, { size: 25, bold: true, color: C.green, align: "center" });
  addText(slide, "Tm", 898, 596, 151, 45, { size: 25, bold: true, color: C.amber, align: "center" });
  addText(slide, "No binding, epitope-competition or thermal improvement is claimed before measurement.", 66, 654, 1145, 26, { size: 17, color: C.muted, valign: "top" });

  slide.speakerNotes.text = [
    "KD — константа диссоциации, измеряемая по SPR; меньший KD означает более сильное связывание. FACS и ridge проверяются на раздельных последовательностных группах. Во внешнем SPR переносе рассмотрено 36 совместимых контролей, интервалы неопределённости не позволяют уверенно выбрать победителя.",
    "ipTM — уверенность модели в интерфейсе, а не измеренная аффинность. PAE показывает неопределённость взаимного положения частей комплекса. Структурный допуск требует одновременно вернуть контакты эпитопа, близко восстановить позу и не иметь тяжёлых столкновений в как минимум двух из трёх seed.",
    "Фактические контроли завершены: WT прошёл 1/3 seed без MSA и 0/3 с target-only MSA при заранее зафиксированном пороге 2/3. Кандидатные refold не запускались, candidate poses не аудировались. Tm — температура плавления белка; термостабильность и конкуренция остаются прямыми лабораторными измерениями.",
    "Sources: docs/research/benchmark-evidence.md; docs/research/external-spr-results.md; docs/research/template-controls-results.md; docs/research/candidate-msa-results.md; docs/research/structure-admission-rule.md; docs/research/target-msa-protocol.md; results/candidate-msa-screen/report.md; docs/research/thermal-assessment.md.",
    "External sources: https://files.rcsb.org/download/1N8Z.cif; https://files.rcsb.org/download/6OGE.cif; https://github.com/AbSciBio/unlocking-de-novo-antibody-design/blob/cd2bd4b6932b5301d435555dc22084d407e4b92b/README.md.",
  ].join("\n\n");
}

if (ppt.slides.items.length !== 2) {
  throw new Error(`Expected exactly 2 slides, got ${ppt.slides.items.length}`);
}

const expectedSlideSizeEmu = "12192000,6858000";
await (await PresentationFile.exportPptx(ppt)).save(candidatePath);

const requirements = {
  explicitTotalSlideCount: 2,
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
if (finalDeck.slides.items.length !== 2) {
  throw new Error(`Final PPTX has ${finalDeck.slides.items.length} slides, expected 2.`);
}
for (const [index, slide] of finalDeck.slides.items.entries()) {
  const image = await finalDeck.export({ slide, format: "png", scale: 1 });
  await fs.writeFile(path.join(outputDir, `slide-${String(index + 1).padStart(2, "0")}.png`), new Uint8Array(await image.arrayBuffer()));
}

console.log(JSON.stringify({
  slideCount: finalDeck.slides.items.length,
  finalPath,
  previews: [path.join(outputDir, "slide-01.png"), path.join(outputDir, "slide-02.png")],
  validationReceipt: receiptPath,
}, null, 2));
