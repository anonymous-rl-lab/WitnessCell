import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const experimentDir = path.dirname(fileURLToPath(import.meta.url));
const resultDir = path.join(experimentDir, "results/formal_30split");
const outputPath = path.join(resultDir, "WitnessCell_Gate07_Norman_results.xlsx");
const previewDir = process.env.WITNESSCELL_PREVIEW_DIR || "/tmp/witnesscell_gate07_previews";
await fs.mkdir(previewDir, { recursive: true });

const workbook = Workbook.create();
const summary = workbook.worksheets.add("Summary");

function parseCsv(text) {
  return text.trimEnd().split(/\r?\n/).map((line) => line.split(",").map((cell) => {
    if (cell === "") return null;
    const numeric = Number(cell);
    return Number.isFinite(numeric) ? numeric : cell;
  }));
}

for (const [file, sheetName] of [
  ["per_seed_metrics.csv", "PerSeed"],
  ["risk_alignment.csv", "RiskAlignment"],
  ["confidence_intervals.csv", "ConfidenceIntervals"],
  ["paired_tests.csv", "PairedTests"],
]) {
  const csv = await fs.readFile(path.join(resultDir, file), "utf8");
  const values = parseCsv(csv);
  const sheet = workbook.worksheets.add(sheetName);
  sheet.getRangeByIndexes(0, 0, values.length, values[0].length).values = values;
}

summary.showGridLines = false;
summary.getRange("A1:F1").merge();
summary.getRange("A1").values = [["WitnessCell Gate 07 — Estimated Witness Risk on Norman"]];
summary.getRange("A1:F1").format = {
  fill: "#17365D",
  font: { bold: true, color: "#FFFFFF", size: 16 },
  horizontalAlignment: "center",
  verticalAlignment: "center",
};
summary.getRange("A1:F1").format.rowHeight = 30;

summary.getRange("A3:C3").values = [["Prediction", "Estimated Witness", "Geometry Only"]];
summary.getRange("A4:A7").values = [
  ["Residual MSE"],
  ["Residual cosine"],
  ["Full-effect cosine"],
  ["Selected discrepancy weight (rho)"],
];
summary.getRange("B4:B7").formulas = [
  ["=AVERAGEIF('PerSeed'!$B$2:$B$91,\"estimated_witness\",'PerSeed'!$C$2:$C$91)"],
  ["=AVERAGEIF('PerSeed'!$B$2:$B$91,\"estimated_witness\",'PerSeed'!$E$2:$E$91)"],
  ["=AVERAGEIF('PerSeed'!$B$2:$B$91,\"estimated_witness\",'PerSeed'!$H$2:$H$91)"],
  ["=AVERAGEIF('PerSeed'!$B$2:$B$91,\"estimated_witness\",'PerSeed'!$J$2:$J$91)"],
];
summary.getRange("C4:C7").formulas = [
  ["=AVERAGEIF('PerSeed'!$B$2:$B$91,\"geometry_only\",'PerSeed'!$C$2:$C$91)"],
  ["=AVERAGEIF('PerSeed'!$B$2:$B$91,\"geometry_only\",'PerSeed'!$E$2:$E$91)"],
  ["=AVERAGEIF('PerSeed'!$B$2:$B$91,\"geometry_only\",'PerSeed'!$H$2:$H$91)"],
  ["=AVERAGEIF('PerSeed'!$B$2:$B$91,\"geometry_only\",'PerSeed'!$J$2:$J$91)"],
];

summary.getRange("E3:F3").values = [["Primary Gate", "Value"]];
summary.getRange("E4:E9").values = [
  ["Residual MSE reduction"],
  ["95% CI low"],
  ["95% CI high"],
  ["Risk vs oracle Spearman"],
  ["Risk vs realized MSE Spearman"],
  ["Formal verdict"],
];
summary.getRange("F4").formulas = [["='ConfidenceIntervals'!C2"]];
summary.getRange("F5:F8").formulas = [
  ["='ConfidenceIntervals'!E2"],
  ["='ConfidenceIntervals'!F2"],
  ["=AVERAGE('RiskAlignment'!B2:B31)"],
  ["=AVERAGE('RiskAlignment'!C2:C31)"],
];
summary.getRange("F9").values = [["PASS"]];

summary.getRange("A11:F11").merge();
summary.getRange("A11").values = [["Frozen protocol: smoke seeds 0–4; untouched formal seeds 100–129; outer-test double outcomes excluded from K, k_t and k_tt estimation."]];
summary.getRange("A11:F11").format = {
  fill: "#EAF2F8",
  font: { italic: true, color: "#1F4E78" },
  wrapText: true,
};
summary.getRange("A11:F11").format.rowHeight = 42;

for (const range of ["A3:C3", "E3:F3"]) {
  summary.getRange(range).format = {
    fill: "#5B9BD5",
    font: { bold: true, color: "#FFFFFF" },
    horizontalAlignment: "center",
    borders: { preset: "outside", style: "thin", color: "#9EADBA" },
  };
}
summary.getRange("A4:C7").format.borders = { preset: "inside", style: "thin", color: "#D9E2F3" };
summary.getRange("E4:F9").format.borders = { preset: "inside", style: "thin", color: "#D9E2F3" };
summary.getRange("B4:C4").format.numberFormat = "0.0000";
summary.getRange("B5:C7").format.numberFormat = "0.000";
summary.getRange("F4:F8").format.numberFormat = "0.0%";
summary.getRange("F7:F8").format.numberFormat = "0.000";
summary.getRange("F9").format = {
  fill: "#C6E0B4",
  font: { bold: true, color: "#375623" },
  horizontalAlignment: "center",
};
summary.getRange("A1:F11").format.font.name = "Aptos";
summary.getRange("A1:F11").format.verticalAlignment = "center";
summary.getRange("A:A").format.columnWidth = 31;
summary.getRange("B:C").format.columnWidth = 18;
summary.getRange("D:D").format.columnWidth = 4;
summary.getRange("E:E").format.columnWidth = 31;
summary.getRange("F:F").format.columnWidth = 18;

for (const sheetName of ["PerSeed", "RiskAlignment", "ConfidenceIntervals", "PairedTests"]) {
  const sheet = workbook.worksheets.getItem(sheetName);
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(1);
  const used = sheet.getUsedRange();
  used.format.font.name = "Aptos";
  used.format.font.size = 10;
  used.format.verticalAlignment = "center";
  used.format.borders = { preset: "inside", style: "thin", color: "#E7E6E6" };
  used.getRow(0).format = {
    fill: "#17365D",
    font: { bold: true, color: "#FFFFFF", name: "Aptos", size: 10 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
    wrapText: true,
  };
  used.getRow(0).format.rowHeight = 32;
  used.format.autofitColumns();
  used.format.autofitRows();
}

workbook.worksheets.getItem("PerSeed").getRange("C2:N91").format.numberFormat = "0.0000";
workbook.worksheets.getItem("RiskAlignment").getRange("B2:G31").format.numberFormat = "0.0000";
workbook.worksheets.getItem("RiskAlignment").getRange("E2:E31").format.numberFormat = "0.00E+00";
workbook.worksheets.getItem("ConfidenceIntervals").getRange("C2:F8").format.numberFormat = "0.0000";
workbook.worksheets.getItem("PairedTests").getRange("C2:G7").format.numberFormat = "0.0000";
workbook.worksheets.getItem("PairedTests").getRange("G2:G7").format.numberFormat = "0.00E+00";
workbook.worksheets.getItem("PairedTests").getRange("D:D").format.columnWidth = 25;
workbook.worksheets.getItem("PairedTests").getRange("G:G").format.columnWidth = 22;

for (const [sheetName, range] of [
  ["Summary", "A1:F11"],
  ["PerSeed", "A1:N20"],
  ["RiskAlignment", "A1:G20"],
  ["ConfidenceIntervals", "A1:G8"],
  ["PairedTests", "A1:G7"],
]) {
  const preview = await workbook.render({ sheetName, range, scale: 1.25, format: "png" });
  await fs.writeFile(
    path.join(previewDir, `${sheetName}.png`),
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const inspection = await workbook.inspect({
  kind: "table",
  sheetId: "Summary",
  range: "A1:F11",
  include: "values,formulas",
  tableMaxRows: 12,
  tableMaxCols: 8,
});
console.log(inspection.ndjson);
const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "final formula error scan",
});
console.log(errors.ndjson);

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(outputPath);
