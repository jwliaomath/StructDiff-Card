import * as $3Dmol from "3dmol";
import "./styles.css";

const data = await fetch("./data/result.json").then((response) => {
  if (!response.ok) throw new Error(`Could not load example result: ${response.status}`);
  return response.json();
});

const residualStops = [
  [0, "#2166ac"],
  [1, "#67a9cf"],
  [2, "#d1e5f0"],
  [3, "#fddbc7"],
  [4, "#ef8a62"],
  [5, "#b2182b"],
];

function residualColor(value) {
  const bounded = Math.max(0, Math.min(5, value));
  const lowerIndex = Math.min(4, Math.floor(bounded));
  const fraction = bounded - lowerIndex;
  const channels = residualStops.slice(lowerIndex, lowerIndex + 2).map(([, hex]) =>
    [1, 3, 5].map((start) => Number.parseInt(hex.slice(start, start + 2), 16)),
  );
  const mixed = channels[0].map((channel, index) =>
    Math.round(channel + (channels[1][index] - channel) * fraction),
  );
  return `#${mixed.map((channel) => channel.toString(16).padStart(2, "0")).join("")}`;
}

function viewerResidueId(value) {
  return /^-?\d+$/.test(value) ? Number.parseInt(value, 10) : value;
}

function fillMetrics() {
  const summary = data.summary;
  document.querySelector("#metric-tm").textContent = summary.tm_reference.toFixed(3);
  document.querySelector("#metric-rmsd").innerHTML = `${summary.rmsd.toFixed(2)} <small>Å</small>`;
  document.querySelector("#metric-aligned").innerHTML = `${summary.aligned_length} <small>/ ${summary.reference_length}</small>`;
  document.querySelector("#metric-chains").textContent = data.chain_alignments.length;
}

function fillChainTable() {
  const body = document.querySelector("#chain-table");
  body.innerHTML = data.chain_alignments
    .map(
      (row) => `<tr>
        <td><span class="chain mobile-chain">${row.mobile_chain}</span></td>
        <td><span class="mapping-arrow">→</span><span class="chain reference-chain">${row.reference_chain}</span></td>
        <td>${row.tm_reference.toFixed(3)}</td>
        <td>${row.rmsd.toFixed(2)} Å</td>
      </tr>`,
    )
    .join("");
  const unmatched = document.querySelector("#unmatched-callout");
  unmatched.innerHTML = `<span>${data.unmatched_mobile_chains.length}</span><div><strong>Unmatched mobile chains</strong><p>${data.unmatched_mobile_chains.join(", ") || "None"}</p></div>`;
}

function fillDifferenceSummary() {
  const summary = data.difference_summary;
  if (!summary) return;
  document.querySelector("#difference-median").textContent = `${summary.median_displacement.toFixed(2)} Å`;
  document.querySelector("#difference-p90").textContent = `${summary.p90_displacement.toFixed(2)} Å`;
  document.querySelector("#difference-max").textContent = `${summary.max_displacement.toFixed(2)} Å`;
  document.querySelector("#difference-threshold-label").textContent = `Residues ≥ ${summary.hotspot_threshold.toFixed(2)} Å`;
  document.querySelector("#difference-high").textContent = `${summary.residues_above_threshold} · ${(summary.fraction_above_threshold * 100).toFixed(1)}%`;
  document.querySelector("#hotspot-table").innerHTML = data.difference_regions
    .map((region) => `<tr>
      <td><strong>${region.rank}</strong></td>
      <td>${region.mobile_chain} → ${region.reference_chain}</td>
      <td>${region.mobile_start}–${region.mobile_end}</td>
      <td>${region.reference_start}–${region.reference_end}</td>
      <td>${region.mean_displacement.toFixed(2)} Å</td>
      <td><strong>${region.peak_displacement.toFixed(2)} Å</strong></td>
    </tr>`)
    .join("");
}

async function setupViewer() {
  const element = document.querySelector("#molecule-viewer");
  const viewer = $3Dmol.createViewer(element, { backgroundColor: "#07131b", antialias: true });
  const [referencePdb, mobilePdb] = await Promise.all([
    fetch("./structures/reference.pdb").then((response) => response.text()),
    fetch("./structures/mobile_aligned_to_reference.pdb").then((response) => response.text()),
  ]);
  viewer.addModel(referencePdb, "pdb");
  viewer.addModel(mobilePdb, "pdb");

  function styleViewer(showUnmatched = true, showDifference = true) {
    viewer.setStyle({}, {});
    viewer.removeAllShapes();
    viewer.removeAllLabels();
    viewer.setStyle(
      { model: 0 },
      { cartoon: { color: showDifference ? "#d8e2e6" : "#ff8a5c", opacity: showDifference ? 0.30 : 0.88 } },
    );
    viewer.setStyle({ model: 0, resn: "HEM" }, { stick: { colorscheme: "Jmol", radius: 0.2 } });
    ["C", "D"].forEach((chain) => {
      viewer.setStyle({ model: 1, chain }, { cartoon: { color: "#42d5c5", opacity: 0.78 } });
    });
    if (showDifference) {
      data.residue_mapping
        .filter((row) => row.displacement !== null)
        .forEach((row) => {
          viewer.setStyle(
            { model: 1, chain: row.mobile_chain, resi: viewerResidueId(row.mobile_residue) },
            { cartoon: { color: residualColor(row.displacement), opacity: 0.92 } },
          );
        });
      data.spatial_patches
        .filter((patch) => patch.member_count >= 2)
        .slice(0, 8)
        .forEach((patch) => {
          const center = { x: patch.centroid_x, y: patch.centroid_y, z: patch.centroid_z };
          viewer.addSphere({ center, radius: 0.65, color: "#ffcf56", opacity: 0.82 });
          viewer.addLabel(`P${patch.rank}`, {
            position: center,
            backgroundColor: "#07131b",
            backgroundOpacity: 0.78,
            fontColor: "#ffcf56",
            fontSize: 11,
            inFront: true,
          });
        });
      const threshold = data.difference_summary?.hotspot_threshold ?? 2;
      data.residue_mapping
        .filter((row) => row.displacement !== null && row.displacement >= threshold)
        .sort((left, right) => right.displacement - left.displacement)
        .slice(0, 20)
        .forEach((row) => {
          const mobileSelection = {
            model: 1,
            chain: row.mobile_chain,
            resi: viewerResidueId(row.mobile_residue),
          };
          const referenceSelection = {
            model: 0,
            chain: row.reference_chain,
            resi: viewerResidueId(row.reference_residue),
          };
          const mobileAtoms = viewer.getModel(1).selectedAtoms(mobileSelection);
          const referenceAtoms = viewer.getModel(0).selectedAtoms(referenceSelection);
          const mobileAnchor = mobileAtoms.find((atom) => atom.atom === "CA" || atom.atom === "C3'");
          const referenceAnchor = referenceAtoms.find((atom) => atom.atom === "CA" || atom.atom === "C3'");
          if (mobileAnchor && referenceAnchor) {
            viewer.addLine({
              start: { x: referenceAnchor.x, y: referenceAnchor.y, z: referenceAnchor.z },
              end: { x: mobileAnchor.x, y: mobileAnchor.y, z: mobileAnchor.z },
              color: "#ffcf56",
              dashed: true,
            });
          }
        });
    }
    if (showUnmatched) {
      ["A", "B"].forEach((chain) => {
        viewer.setStyle({ model: 1, chain }, { cartoon: { color: "#b77ce8", opacity: 0.4 } });
      });
    }
    viewer.setStyle({ model: 1, resn: "HEM" }, { stick: { colorscheme: "Jmol", radius: 0.17 } });
    viewer.render();
  }

  styleViewer(true);
  viewer.zoomTo();
  viewer.zoom(1.05, 400);

  const unmatchedToggle = document.querySelector("#toggle-unmatched");
  const differenceToggle = document.querySelector("#toggle-difference");
  const updateStyle = () => styleViewer(unmatchedToggle.checked, differenceToggle.checked);
  unmatchedToggle.addEventListener("change", updateStyle);
  differenceToggle.addEventListener("change", updateStyle);
  document.querySelector("#reset-view").addEventListener("click", () => {
    viewer.zoomTo();
    viewer.zoom(1.05, 300);
  });
  const spinButton = document.querySelector("#spin-button");
  let spinning = false;
  spinButton.addEventListener("click", () => {
    spinning = !spinning;
    viewer.spin(spinning ? "y" : false);
    spinButton.textContent = spinning ? "Stop" : "Spin";
  });
  return viewer;
}

function drawDisplacementChart() {
  const svg = document.querySelector("#displacement-chart");
  const tooltip = document.querySelector("#chart-tooltip");
  const rows = data.residue_mapping.filter((row) => row.displacement !== null);
  const width = 1120;
  const height = 390;
  const margin = { top: 28, right: 24, bottom: 54, left: 66 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const maxValue = Math.max(5.5, ...rows.map((row) => row.displacement)) * 1.08;
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);

  const x = (index) => margin.left + (index / Math.max(1, rows.length - 1)) * plotWidth;
  const y = (value) => margin.top + plotHeight - (value / maxValue) * plotHeight;
  const ns = "http://www.w3.org/2000/svg";
  const add = (name, attrs = {}, text = "") => {
    const node = document.createElementNS(ns, name);
    Object.entries(attrs).forEach(([key, value]) => node.setAttribute(key, value));
    if (text) node.textContent = text;
    svg.appendChild(node);
    return node;
  };

  [0, 1, 2, 3, 4, 5].forEach((value) => {
    add("line", { x1: margin.left, x2: width - margin.right, y1: y(value), y2: y(value), class: "guide" });
    add("text", { x: margin.left - 13, y: y(value) + 5, class: "axis-label", "text-anchor": "end" }, `${value}`);
  });
  const threshold = data.difference_summary?.hotspot_threshold ?? 2;
  add("line", { x1: margin.left, x2: width - margin.right, y1: y(threshold), y2: y(threshold), class: "guide threshold" });
  add("text", { x: width - margin.right - 5, y: y(threshold) - 8, class: "threshold-label", "text-anchor": "end" }, `hotspot threshold ${threshold.toFixed(2)} Å`);
  add("text", { x: 18, y: height / 2, class: "axis-title", transform: `rotate(-90 18 ${height / 2})`, "text-anchor": "middle" }, "Displacement (Å)");
  add("text", { x: width / 2, y: height - 12, class: "axis-title", "text-anchor": "middle" }, "Matched residue pairs in chain-mapping order");

  const grouped = new Map();
  rows.forEach((row, index) => {
    const key = `${row.mobile_chain}→${row.reference_chain}`;
    if (!grouped.has(key)) grouped.set(key, []);
    grouped.get(key).push({ row, index });
  });
  const colors = ["#0d9488", "#287fb8", "#b45fcc", "#d76542"];
  [...grouped.entries()].forEach(([key, values], groupIndex) => {
    add("polyline", {
      points: values.map(({ row, index }) => `${x(index)},${y(row.displacement)}`).join(" "),
      fill: "none",
      stroke: "#b9c8c6",
      "stroke-width": 1.4,
      "stroke-linejoin": "round",
      class: "data-line",
    });
    const radius = 3;
    values.forEach(({ row, index }) => {
      add("circle", { cx: x(index), cy: y(row.displacement), r: radius, fill: residualColor(row.displacement), class: "data-point" });
    });
    const smoothed = values.map(({ row }, valueIndex) => {
      const start = Math.max(0, valueIndex - 3);
      const end = Math.min(values.length, valueIndex + 4);
      return values.slice(start, end).reduce((sum, item) => sum + item.row.displacement, 0) / (end - start);
    });
    add("polyline", {
      points: values.map(({ index }, valueIndex) => `${x(index)},${y(smoothed[valueIndex])}`).join(" "),
      fill: "none",
      stroke: colors[groupIndex % colors.length],
      "stroke-width": 3,
      "stroke-linejoin": "round",
      class: "smooth-line",
    });
    const first = values[0];
    add("text", { x: x(first.index), y: margin.top + 14 + groupIndex * 20, fill: colors[groupIndex % colors.length], class: "series-label" }, key);
  });

  rows.forEach((row, index) => {
    const hit = add("circle", { cx: x(index), cy: y(row.displacement), r: 7, class: "chart-hit" });
    hit.addEventListener("mouseenter", (event) => {
      tooltip.hidden = false;
      tooltip.innerHTML = `<strong>${row.mobile_chain}:${row.mobile_residue} ↔ ${row.reference_chain}:${row.reference_residue}</strong><span>${row.displacement.toFixed(2)} Å</span>`;
      const bounds = svg.getBoundingClientRect();
      tooltip.style.left = `${event.clientX - bounds.left + 14}px`;
      tooltip.style.top = `${event.clientY - bounds.top - 48}px`;
    });
    hit.addEventListener("mouseleave", () => { tooltip.hidden = true; });
  });
  data.difference_regions.forEach((region) => {
    const peakIndex = rows.findIndex(
      (row) => row.mobile_chain === region.mobile_chain && row.mobile_residue === region.peak_mobile_residue,
    );
    if (peakIndex < 0) return;
    add("circle", { cx: x(peakIndex), cy: y(region.peak_displacement), r: 6, class: "hotspot-peak" });
    add("text", { x: x(peakIndex), y: y(region.peak_displacement) - 12, class: "hotspot-label", "text-anchor": "middle" }, `${region.rank}`);
  });
}

function fillGapTrack() {
  const root = document.querySelector("#gap-track");
  const rows = data.chain_alignments.map((alignment) => {
    const cells = [...alignment.mobile_sequence].map((mobile, index) => {
      const reference = alignment.reference_sequence[index];
      let status = "matched";
      let label = mobile;
      if (mobile === "-") { status = "only-reference"; label = reference; }
      else if (reference === "-") status = "only-mobile";
      return `<i class="${status}" title="${status.replace("-", " ")}">${label}</i>`;
    }).join("");
    return `<div class="sequence-row"><div><strong>${alignment.mobile_chain} → ${alignment.reference_chain}</strong><span>${alignment.aligned_length} aligned · ${(alignment.coverage_reference * 100).toFixed(0)}% reference coverage</span></div><div class="sequence-cells">${cells}</div></div>`;
  }).join("");
  const unmatched = data.unmatched_mobile_chains.map((chain) => `<span>Mobile chain ${chain}<b>unmatched</b></span>`).join("");
  root.innerHTML = `${rows}<div class="unmatched-track"><strong>Unmatched-chain track</strong><div>${unmatched || "None"}</div></div>`;
}

function drawContactMap() {
  const summary = data.contact_summary;
  document.querySelector("#contact-gained").textContent = summary.gained_contacts;
  document.querySelector("#contact-lost").textContent = summary.lost_contacts;
  document.querySelector("#contact-interchain").textContent =
    summary.inter_chain_gained_contacts + summary.inter_chain_lost_contacts;

  const svg = document.querySelector("#contact-map");
  const width = 560;
  const height = 400;
  const margin = 48;
  const maximum = Math.max(
    1,
    ...data.residue_mapping.map((row) => row.alignment_index),
  );
  const x = (value) => margin + (value / maximum) * (width - margin * 1.4);
  const y = (value) => height - margin - (value / maximum) * (height - margin * 1.4);
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.innerHTML = `
    <line x1="${margin}" y1="${height - margin}" x2="${width - 20}" y2="${height - margin}" class="contact-axis" />
    <line x1="${margin}" y1="${height - margin}" x2="${margin}" y2="20" class="contact-axis" />
    <line x1="${margin}" y1="${height - margin}" x2="${width - 20}" y2="20" class="contact-diagonal" />
    <text x="${width / 2}" y="390" class="contact-axis-label" text-anchor="middle">mapped residue index 1</text>
    <text x="14" y="${height / 2}" class="contact-axis-label" text-anchor="middle" transform="rotate(-90 14 ${height / 2})">mapped residue index 2</text>`;
  const namespace = "http://www.w3.org/2000/svg";
  data.contact_changes.forEach((row) => {
    const point = document.createElementNS(namespace, row.change === "gained" ? "circle" : "path");
    const pointX = x(row.alignment_index_1);
    const pointY = y(row.alignment_index_2);
    if (row.change === "gained") {
      point.setAttribute("cx", pointX);
      point.setAttribute("cy", pointY);
      point.setAttribute("r", row.inter_chain ? 5 : 3.5);
      point.setAttribute("class", "contact-point gained");
    } else {
      const radius = row.inter_chain ? 5 : 3.5;
      point.setAttribute(
        "d",
        `M ${pointX - radius} ${pointY - radius} L ${pointX + radius} ${pointY + radius} M ${pointX + radius} ${pointY - radius} L ${pointX - radius} ${pointY + radius}`,
      );
      point.setAttribute("class", "contact-point lost");
    }
    const title = document.createElementNS(namespace, "title");
    title.textContent = `${row.change}: ${row.mobile_chain_1}:${row.mobile_residue_1} — ${row.mobile_chain_2}:${row.mobile_residue_2}; mobile ${row.mobile_distance.toFixed(2)} Å, reference ${row.reference_distance.toFixed(2)} Å`;
    point.appendChild(title);
    svg.appendChild(point);
  });
  if (!data.contact_changes.length) {
    const message = document.createElementNS(namespace, "text");
    message.setAttribute("x", width / 2);
    message.setAttribute("y", height / 2);
    message.setAttribute("class", "contact-empty");
    message.setAttribute("text-anchor", "middle");
    message.textContent = "No stable contact changes detected";
    svg.appendChild(message);
  }
}

function fillMotionDecomposition() {
  const maximum = Math.max(
    1,
    ...data.motion_decomposition.flatMap((row) => [row.rigid_body_rmsd, row.internal_rmsd]),
  );
  document.querySelector("#motion-bars").innerHTML = data.motion_decomposition
    .map((row) => `<div class="motion-row">
      <strong>${row.mobile_chain}→${row.reference_chain}</strong>
      <div><i class="rigid" style="width:${(row.rigid_body_rmsd / maximum) * 100}%"></i><span>${row.rigid_body_rmsd.toFixed(2)} Å rigid</span></div>
      <div><i class="internal" style="width:${(row.internal_rmsd / maximum) * 100}%"></i><span>${row.internal_rmsd.toFixed(2)} Å internal</span></div>
    </div>`)
    .join("");
  document.querySelector("#motion-table").innerHTML = data.motion_decomposition
    .map((row) => `<tr>
      <td><strong>${row.mobile_chain}→${row.reference_chain}</strong></td>
      <td>${row.rigid_body_rmsd.toFixed(2)} Å</td>
      <td>${row.internal_rmsd.toFixed(2)} Å</td>
      <td>${row.classification}</td>
    </tr>`)
    .join("");
}

function fillSpatialPatches() {
  document.querySelector("#patch-table").innerHTML = data.spatial_patches
    .map((patch) => `<tr title="${patch.mobile_residues}">
      <td><strong>P${patch.rank}</strong></td>
      <td>${patch.member_count}</td>
      <td>${patch.mobile_chains} → ${patch.reference_chains}</td>
      <td>${patch.mean_displacement.toFixed(2)} Å</td>
      <td><strong>${patch.peak_displacement.toFixed(2)} Å</strong></td>
      <td>${patch.spatial_radius.toFixed(2)} Å</td>
    </tr>`)
    .join("");
}

function setupCopyButton() {
  const button = document.querySelector("#copy-command");
  button.addEventListener("click", async () => {
    const command = document.querySelector("#docker-command").textContent;
    await navigator.clipboard.writeText(command);
    button.textContent = "Copied";
    setTimeout(() => { button.textContent = "Copy"; }, 1500);
  });
}

fillMetrics();
fillChainTable();
fillDifferenceSummary();
drawDisplacementChart();
drawContactMap();
fillMotionDecomposition();
fillSpatialPatches();
fillGapTrack();
setupCopyButton();
setupViewer().catch((error) => {
  document.querySelector("#molecule-viewer").innerHTML = `<p class="viewer-error">3D viewer could not load: ${error.message}</p>`;
});
