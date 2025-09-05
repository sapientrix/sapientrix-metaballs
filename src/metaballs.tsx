import React, { useEffect, useMemo, useRef, useState } from "react";

// 3x3 Metaballs Grid using the SVG membrane approach (no filters, no pixels)
// Controls: radius, gap, min connection thickness (v), stroke thickness, colors, download SVG.
// Click behavior: clicking a CIRCLE or a METABALL toggles 3 states:
//   0 → path = pathColor (stroke only)
//   1 → path = backgroundColor (effectively hidden)
//   2 → filled = pathColor (fill the shape; outlines may still render depending on other toggles)

// ---------------- Constants ----------------
const XML_HEADER = '<?xml version="1.0" encoding="UTF-8"?>\n'; // used by downloadSVG

// ---------------- Utils (geometry) ----------------
const TAU = Math.PI * 2;
const EPS = 1e-6;

const dist = (a: [number, number], b: [number, number]) => {
  const dx = b[0] - a[0];
  const dy = b[1] - a[1];
  return Math.hypot(dx, dy);
};

const angle = (to: [number, number], from: [number, number]) => {
  const dx = to[0] - from[0];
  const dy = to[1] - from[1];
  return Math.atan2(dy, dx);
};

const getVector = (
  origin: [number, number],
  ang: number,
  len: number
): [number, number] => {
  return [origin[0] + Math.cos(ang) * len, origin[1] + Math.sin(ang) * len];
};

const pt = (p: [number, number]) => `${p[0].toFixed(3)},${p[1].toFixed(3)}`;
const norm = (a: number) => {
  let t = a % TAU;
  if (t < 0) t += TAU;
  return t;
};

// cubic Bezier evaluator and sampler
function cubicAt(
  p0: number[],
  p1: number[],
  p2: number[],
  p3: number[],
  t: number
): [number, number] {
  const mt = 1 - t;
  const a = mt * mt * mt;
  const b = 3 * mt * mt * t;
  const c = 3 * mt * t * t;
  const d = t * t * t;
  return [
    a * p0[0] + b * p1[0] + c * p2[0] + d * p3[0],
    a * p0[1] + b * p1[1] + c * p2[1] + d * p3[1],
  ];
}

function sampleCubic(
  p0: [number, number],
  c1: [number, number],
  c2: [number, number],
  p1: [number, number],
  n = 48
) {
  const pts: [number, number][] = [];
  for (let i = 0; i <= n; i++) {
    const t = i / n;
    pts.push(cubicAt(p0, c1, c2, p1, t));
  }
  return pts;
}

// Return the minor (<= PI) CCW arc span from a to b as [start,end] with end >= start.
function minorSpan(a: number, b: number): [number, number] {
  let s = norm(a);
  let e = norm(b);
  let d = (e - s + TAU) % TAU; // CCW distance from s to e
  if (d > Math.PI) {
    // take the other way around so it is <= PI
    s = e;
    d = (norm(a) - norm(b) + TAU) % TAU; // now this is <= PI
  }
  return [s, s + d];
}

// Merge overlapping [start,end] spans (both within [0, 2π] and end >= start)
function mergeSpans(spans: Array<[number, number]>): Array<[number, number]> {
  if (spans.length === 0) return [];
  const sorted = [...spans].sort((a, b) => a[0] - b[0]);
  const merged: Array<[number, number]> = [];
  let [cs, ce] = sorted[0];
  for (let i = 1; i < sorted.length; i++) {
    const [s, e] = sorted[i];
    if (s <= ce + EPS) {
      ce = Math.max(ce, e);
    } else {
      merged.push([cs, ce]);
      [cs, ce] = [s, e];
    }
  }
  merged.push([cs, ce]);
  return merged;
}

// Compute complement segments of [0,2π] given minor spans to remove
function complementSegments(
  spans: Array<[number, number]>
): Array<[number, number]> {
  if (spans.length === 0) return [[0, TAU]];
  // Split wraps
  const parts: Array<[number, number]> = [];
  for (const [s0, e0] of spans) {
    const s = s0;
    const e = e0;
    if (e <= TAU) {
      parts.push([s, e]);
    } else {
      parts.push([s, TAU]);
      parts.push([0, e - TAU]);
    }
  }
  const merged = mergeSpans(parts);
  const complement: Array<[number, number]> = [];
  let prev = 0;
  for (const [s, e] of merged) {
    if (s > prev + EPS) complement.push([prev, s]);
    prev = e;
  }
  if (prev < TAU - EPS) complement.push([prev, TAU]);
  return complement;
}

// Build an SVG arc command from start→end CCW (sweep=1) on a circle
function arcCmd(
  center: [number, number],
  r: number,
  start: number,
  end: number
): string {
  const pStart = getVector(center, start, r);
  const pEnd = getVector(center, end, r);
  const delta = end - start;
  const largeArc = delta > Math.PI ? 1 : 0;
  const sweep = 1; // CCW in our angle convention
  return [
    "M",
    pt(pStart),
    "A",
    r.toFixed(3),
    r.toFixed(3),
    "0",
    String(largeArc),
    String(sweep),
    pt(pEnd),
  ].join(" ");
}

function arcPoints(
  center: [number, number],
  r: number,
  start: number,
  end: number,
  n = 32
) {
  const [s, e] = minorSpan(start, end);
  const span = e - s;
  const steps = Math.max(8, Math.ceil((span * r) / 8));
  const pts: [number, number][] = [];
  for (let i = 1; i <= steps; i++) {
    const a = s + (span * i) / steps;
    pts.push(getVector(center, a, r));
  }
  return pts;
}

// Simple point-in-polygon (even-odd)
function pointInPolygon(p: [number, number], poly: [number, number][]) {
  let inside = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const xi = poly[i][0],
      yi = poly[i][1];
    const xj = poly[j][0],
      yj = poly[j][1];
    const intersect =
      yi > p[1] !== yj > p[1] &&
      p[0] < ((xj - xi) * (p[1] - yi)) / (yj - yi + 1e-12) + xi;
    if (intersect) inside = !inside;
  }
  return inside;
}

const pathFromPoints = (pts: [number, number][]) =>
  pts.length
    ? `M ${pt(pts[0])} ` +
      pts
        .slice(1)
        .map((q) => `L ${pt(q)}`)
        .join(" ")
    : "";
const polygonPath = (pts: [number, number][]) =>
  pts.length
    ? `M ${pt(pts[0])} ` +
      pts
        .slice(1)
        .map((q) => `L ${pt(q)}`)
        .join(" ") +
      " Z"
    : "";

// Distance from point to segment
function segDist(p: [number, number], a: [number, number], b: [number, number]) {
  const vx = b[0] - a[0];
  const vy = b[1] - a[1];
  const wx = p[0] - a[0];
  const wy = p[1] - a[1];
  const len2 = vx * vx + vy * vy + 1e-12;
  const t = Math.max(0, Math.min(1, (wx * vx + wy * vy) / len2));
  const cx = a[0] + vx * t;
  const cy = a[1] + vy * t;
  const dx = p[0] - cx;
  const dy = p[1] - cy;
  return Math.hypot(dx, dy);
}

function polylineDist(p: [number, number], poly: [number, number][], closed = true) {
  if (poly.length < 2) return Infinity;
  let best = Infinity;
  for (let i = 0; i < poly.length - 1; i++) {
    best = Math.min(best, segDist(p, poly[i], poly[i + 1]));
  }
  if (closed) best = Math.min(best, segDist(p, poly[poly.length - 1], poly[0]));
  return best;
}

/**
 * Membrane parameters between two circles using Varun Vachhar's membrane formulation.
 * Returns { path: (two cubic curves), angles for circle1 (a1,a2) and circle2 (a3,a4), and geometry }.
 * If no connection is possible, returns null.
 */
function metaballParams(
  radius1: number,
  radius2: number,
  center1: [number, number],
  center2: [number, number],
  handleSize = 2.4,
  v = 0.5,
  reachFactor = 1.8
): null | {
  path: string;
  a1: number;
  a2: number;
  a3: number;
  a4: number;
  p1: [number, number];
  p2: [number, number];
  p3: [number, number];
  p4: [number, number];
  h1: [number, number];
  h2: [number, number];
  h3: [number, number];
  h4: [number, number];
} {
  const HALF_PI = Math.PI / 2;
  const d = dist(center1, center2);
  // Allow connections up to reachFactor × (r1 + r2)
  const maxDist = (radius1 + radius2) * reachFactor;
  let u1 = 0,
    u2 = 0;

  if (
    radius1 === 0 ||
    radius2 === 0 ||
    d > maxDist ||
    d <= Math.abs(radius1 - radius2)
  ) {
    return null;
  }

  if (d < radius1 + radius2) {
    u1 = Math.acos(
      (radius1 * radius1 + d * d - radius2 * radius2) / (2 * radius1 * d)
    );
    u2 = Math.acos(
      (radius2 * radius2 + d * d - radius1 * radius1) / (2 * radius2 * d)
    );
  }

  const angBetween = angle(center2, center1);
  const maxSpread = Math.acos((radius1 - radius2) / d);

  const a1 = angBetween + u1 + (maxSpread - u1) * v;
  const a2 = angBetween - u1 - (maxSpread - u1) * v;
  const a3 = angBetween + Math.PI - u2 - (Math.PI - u2 - maxSpread) * v;
  const a4 = angBetween - Math.PI + u2 + (Math.PI - u2 - maxSpread) * v;

  const p1 = getVector(center1, a1, radius1);
  const p2 = getVector(center1, a2, radius1);
  const p3 = getVector(center2, a3, radius2);
  const p4 = getVector(center2, a4, radius2);

  const totalRadius = radius1 + radius2;
  const d2Base = Math.min(v * handleSize, dist(p1, p3) / totalRadius);
  const d2 = d2Base * Math.min(1, (d * 2) / (radius1 + radius2));
  const r1 = radius1 * d2;
  const r2 = radius2 * d2;

  const h1 = getVector(p1, a1 - Math.PI / 2, r1);
  const h2 = getVector(p2, a2 + Math.PI / 2, r1);
  const h3 = getVector(p3, a3 + Math.PI / 2, r2);
  const h4 = getVector(p4, a4 - Math.PI / 2, r2);

  // Two cubic curves (no circle arc) — membrane only
  const path = [
    "M",
    pt(p1),
    "C",
    pt(h1),
    pt(h3),
    pt(p3),
    "M",
    pt(p4),
    "C",
    pt(h4),
    pt(h2),
    pt(p2),
  ].join(" ");

  return { path, a1, a2, a3, a4, p1, p2, p3, p4, h1, h2, h3, h4 };
}

// ---------------- UI ----------------
export default function MetaballsGrid() {
  // URL params helpers
  const readParams = () =>
    typeof window !== "undefined" ? new URLSearchParams(window.location.search) : null;
  const readNum = (name: string, def: number, min?: number, max?: number) => {
    const p = readParams();
    if (!p) return def;
    if (!p.has(name)) return def;
    const raw = p.get(name);
    if (raw == null || raw === "") return def;
    const v = Number(raw);
    if (!Number.isFinite(v)) return def;
    if (min !== undefined && v < min) return def;
    if (max !== undefined && v > max) return def;
    return v;
  };
  const readStr = (name: string, def: string) => {
    const p = readParams();
    if (!p) return def;
    const v = p.get(name);
    return v ? v : def;
  };
  const readBool = (name: string, def: boolean) => {
    const p = readParams();
    if (!p) return def;
    const v = p.get(name);
    if (v == null) return def;
    return v === "1" || v === "true";
  };

  const [radius, setRadius] = useState(readNum("r", 30, 8, 80));
  const [gap, setGap] = useState(readNum("g", 8, 0, 120));
  const [vMin, setVMin] = useState(readNum("v", 0.5, 0, 1));
  const [strokeW, setStrokeW] = useState(readNum("sw", 2, 1, 16));
  const [reach, setReach] = useState(readNum("reach", 1.8, 0.8, 3));
  const [autoReach, setAutoReach] = useState(readBool("ar", true));
  // Advanced tolerances (for trimming + fill overlap)
  const [advOpen, setAdvOpen] = useState(false);
  const autoDelta = () => Math.max(0.75 * strokeW, 1.5);
  const autoOverlap = () => Math.max(0.75 * strokeW, 1.5);
  const autoJoin = () => Math.max(1.25 * strokeW, 3);
  const autoMinFrag = () => Math.max(2.5 * strokeW, 6);
  const [advDelta, setAdvDelta] = useState(readNum("ad", autoDelta(), 0.2, 10));
  const [advOverlap, setAdvOverlap] = useState(readNum("ao", autoOverlap(), 0, 12));
  const [advJoin, setAdvJoin] = useState(readNum("aj", autoJoin(), 0, 30));
  const [advMinFrag, setAdvMinFrag] = useState(readNum("am", autoMinFrag(), 0, 60));
  const [advAuto, setAdvAuto] = useState(readBool("aa", true));
  const [bgColor, setBgColor] = useState(readStr("bc", "#1E1E1E"));
  const [pathColor, setPathColor] = useState(readStr("pc", "#D9D9D9"));

  const svgRef = useRef<SVGSVGElement | null>(null);

  // Grid layout math
  const pad = 40; // outer padding
  const dx = 2 * radius + gap; // center-to-center spacing
  const dy = dx;
  const cols = 3,
    rows = 3;
  const width = pad * 2 + (cols - 1) * dx + 2 * radius;
  const height = pad * 2 + (rows - 1) * dy + 2 * radius;

  // Centers for 3x3 grid (row-major: r=1..3 top→bottom, c=1..3 left→right)
  const centers: Record<string, [number, number]> = useMemo(() => {
    const map: Record<string, [number, number]> = {};
    for (let r = 1; r <= rows; r++) {
      for (let c = 1; c <= cols; c++) {
        const x = pad + (c - 1) * dx + radius;
        const y = pad + (r - 1) * dy + radius;
        map[`${r},${c}`] = [x, y];
      }
    }
    return map;
  }, [dx, dy, pad, radius]);

  // Metaball groups (connect only within each group)
  const group1 = ["1,1", "2,2"]; // (1,1) ↔ (2,2)
  const group2 = ["1,2", "2,3", "3,2"]; // all pairwise connections among these three

  // Precompute connection pairs
  const pairs: Array<[string, string]> = useMemo(() => {
    const p: Array<[string, string]> = [];
    p.push([group1[0], group1[1]]);
    for (let i = 0; i < group2.length; i++) {
      for (let j = i + 1; j < group2.length; j++) {
        p.push([group2[i], group2[j]]);
      }
    }
    return p;
  }, []);

  // Adaptive reach recommendation based on current layout and pairs
  const recommendedReach = useMemo(() => {
    const sumR = radius * 2;
    if (sumR <= 0) return 1.0;
    const diagNorm = Math.hypot(dx, dy) / sumR; // longest single-step
    let maxWithinDiag = 1.0;
    for (const [ka, kb] of pairs) {
      const c1 = centers[ka];
      const c2 = centers[kb];
      if (!c1 || !c2) continue;
      const dn = dist(c1, c2) / sumR;
      if (dn <= diagNorm + 1e-6) maxWithinDiag = Math.max(maxWithinDiag, dn);
    }
    const margin = 0.03; // headroom for numerical stability
    return Math.min(3, Math.max(0.8, maxWithinDiag + margin));
  }, [dx, dy, centers, pairs, radius]);

  useEffect(() => {
    if (autoReach) setReach(recommendedReach);
  }, [autoReach, recommendedReach]);

  useEffect(() => {
    if (advAuto) {
      setAdvDelta(autoDelta());
      setAdvOverlap(autoOverlap());
      setAdvJoin(autoJoin());
      setAdvMinFrag(autoMinFrag());
    }
  }, [advAuto, strokeW]);

  // States per circle and per pair (0=stroke pathColor, 1=stroke bgColor, 2=fill pathColor)
  const [circleStates, setCircleStates] = useState<Record<string, number>>({});
  const [pairStates, setPairStates] = useState<number[]>([]);

  useEffect(() => {
    // initialize when geometry changes
    const next: Record<string, number> = { ...circleStates };
    Object.keys(centers).forEach((k) => {
      if (!(k in next)) next[k] = 0;
    });
    setCircleStates(next);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [centers]);

  useEffect(() => {
    setPairStates((prev) =>
      prev.length === pairs.length ? prev : new Array(pairs.length).fill(0)
    );
  }, [pairs.length]);

  // Compute membrane paths and the arcs to keep for each circle (outer boundary only)
  const { pairArtifacts, arcsByCircle, connectedKeys } = useMemo(() => {
    type PairGeom = {
      keyA: string;
      keyB: string;
      poly: [number, number][]; // closed polygon approximating the capsule/membrane region
      curve1: [number, number][]; // p1->p3 polyline
      curve2: [number, number][]; // p4->p2 polyline
    };

    const interiorSpans: Record<string, Array<[number, number]>> = {};
    const pushSpan = (key: string, a: number, b: number) => {
      const span = minorSpan(a, b);
      if (!interiorSpans[key]) interiorSpans[key] = [];
      interiorSpans[key].push(span);
    };

    const pairGeoms: PairGeom[] = [];

    for (const [ka, kb] of pairs) {
      const c1 = centers[ka];
      const c2 = centers[kb];
      if (!c1 || !c2) continue;
      const mb = metaballParams(radius, radius, c1, c2, 2.4, vMin, reach);
      if (!mb) continue;

      pushSpan(ka, mb.a2, mb.a1);
      pushSpan(kb, mb.a4, mb.a3);

      const curve1 = sampleCubic(mb.p1, mb.h1, mb.h3, mb.p3, 64);
      // Slight overlap on arcs so fills overlap circles (removes hairline seam)
      const OVERLAP = advOverlap;
      const arc2 = arcPoints(c2, radius + OVERLAP, mb.a3, mb.a4, 48);
      const curve2 = sampleCubic(mb.p4, mb.h4, mb.h2, mb.p2, 64);
      const arc1 = arcPoints(c1, radius + OVERLAP, mb.a2, mb.a1, 48);
      const poly = [...curve1, ...arc2, ...curve2, ...arc1];
      pairGeoms.push({ keyA: ka, keyB: kb, poly, curve1, curve2 });
    }

    // Trim curves against all blobs so we only keep outer segments
    const circles = Object.keys(centers).map((k) => ({ c: centers[k], r: radius }));
    // Union membership without tolerance; we'll probe offsets around edges
    const unionContains = (p: [number, number]) => {
      for (let i = 0; i < circles.length; i++) {
        const { c, r } = circles[i];
        if (dist(p, c) <= r + 1e-9) return true; // include boundary
      }
      for (let i = 0; i < pairGeoms.length; i++) {
        const poly = pairGeoms[i].poly;
        if (poly.length && pointInPolygon(p, poly)) return true;
      }
      return false;
    };

    const pairArtifacts = pairGeoms.map((geom, i) => {
      const segPts: [number, number][][] = [];
      const delta = advDelta;
      const isBorderSeg = (a: [number, number], b: [number, number]) => {
        const dx = b[0] - a[0];
        const dy = b[1] - a[1];
        const L = Math.hypot(dx, dy);
        if (L < 1e-6) return false;
        const nx = -dy / L;
        const ny = dx / L;
        const ts = L < 2 * delta ? [0.5] : [0.2, 0.4, 0.6, 0.8];
        let votes = 0;
        for (const t of ts) {
          const mx = a[0] + dx * t;
          const my = a[1] + dy * t;
          const p1: [number, number] = [mx + nx * delta, my + ny * delta];
          const p2: [number, number] = [mx - nx * delta, my - ny * delta];
          const c1 = unionContains(p1);
          const c2 = unionContains(p2);
          if (c1 !== c2) votes++;
        }
        return votes >= Math.ceil(ts.length / 2);
      };

      for (const polyline of [geom.curve1, geom.curve2]) {
        if (polyline.length < 2) continue;
        let segment: [number, number][] = [];
        for (let k = 0; k < polyline.length - 1; k++) {
          const a = polyline[k];
          const b = polyline[k + 1];
          const keep = isBorderSeg(a, b);
          if (keep) {
            if (segment.length === 0) segment.push(a);
            segment.push(b);
          }
          if ((!keep || k === polyline.length - 2) && segment.length > 1) {
            segPts.push(segment);
            segment = [];
          }
        }
      }

      // Merge small gaps between consecutive segments (bridging near joints)
      const JOIN = advJoin;
      const merged: [number, number][][] = [];
      for (const s of segPts) {
        if (merged.length === 0) {
          merged.push(s);
          continue;
        }
        const prev = merged[merged.length - 1];
        const pa = prev[prev.length - 1];
        const pb = s[0];
        if (dist(pa, pb) <= JOIN) {
          // stitch
          if (dist(pa, pb) > 0) prev.push(pb);
          for (let ii = 1; ii < s.length; ii++) prev.push(s[ii]);
        } else {
          merged.push(s);
        }
      }

      const segLen = (pts: [number, number][]) => {
        let L = 0;
        for (let i = 0; i < pts.length - 1; i++) L += dist(pts[i], pts[i + 1]);
        return L;
      };
      const MINLEN = advMinFrag;
      const segs = merged
        .filter((pts) => segLen(pts) >= MINLEN)
        .map((pts) => pathFromPoints(pts));
      return { fillD: polygonPath(geom.poly), segments: segs };
    });

    // Circle arcs (outer boundary only) grouped by circle key
    const arcsByCircle: Record<string, string[]> = {};
    for (const key of Object.keys(interiorSpans)) {
      const center = centers[key];
      const spans = interiorSpans[key];
      const split: Array<[number, number]> = [];
      for (const [s, e] of spans) {
        if (e <= TAU) split.push([s, e]);
        else {
          split.push([s, TAU]);
          split.push([0, e - TAU]);
        }
      }
      const merged = mergeSpans(split);
      const comps = complementSegments(merged);
      for (const [s, e] of comps) {
        if (e - s < EPS) continue;
        if (!arcsByCircle[key]) arcsByCircle[key] = [];
        arcsByCircle[key].push(arcCmd(center, radius, s, e));
      }
    }

    const connectedKeys = new Set(Object.keys(interiorSpans));
    return { pairArtifacts, arcsByCircle, connectedKeys };
  }, [centers, pairs, radius, vMin, reach, strokeW, advDelta, advOverlap, advJoin, advMinFrag]);

  // Persist settings to URL params
  useEffect(() => {
    if (typeof window === "undefined") return;
    const p = new URLSearchParams(window.location.search);
    p.set("r", String(radius));
    p.set("g", String(gap));
    p.set("v", vMin.toFixed(2));
    p.set("sw", String(strokeW));
    p.set("pc", pathColor);
    p.set("bc", bgColor);
    p.set("reach", reach.toFixed(2));
    p.set("ar", autoReach ? "1" : "0");
    p.set("ad", advDelta.toFixed(2));
    p.set("ao", advOverlap.toFixed(2));
    p.set("aj", advJoin.toFixed(2));
    p.set("am", advMinFrag.toFixed(2));
    p.set("aa", advAuto ? "1" : "0");
    const url = window.location.pathname + "?" + p.toString();
    window.history.replaceState(null, "", url);
  }, [radius, gap, vMin, strokeW, pathColor, bgColor, reach, autoReach, advDelta, advOverlap, advJoin, advMinFrag, advAuto]);

  const toggleCircle = (key: string) => {
    setCircleStates((s) => ({ ...s, [key]: ((s[key] ?? 0) + 1) % 3 }));
  };
  const togglePair = (idx: number) => {
    setPairStates((arr) => {
      const next = [...arr];
      next[idx] = ((next[idx] ?? 0) + 1) % 3;
      return next;
    });
  };

  const downloadSVG = () => {
    const svg = svgRef.current;
    if (!svg) return;
    const clone = svg.cloneNode(true) as SVGSVGElement;
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
    clone.setAttribute("width", String(width));
    clone.setAttribute("height", String(height));

    const serializer = new XMLSerializer();
    const source = serializer.serializeToString(clone);
    const blob = new Blob([XML_HEADER + source], {
      type: "image/svg+xml;charset=utf-8",
    });
    const url = URL.createObjectURL(blob);

    const a = document.createElement("a");
    a.href = url;
    a.download = "metaballs-grid.svg";
    document.body.appendChild(a);
    a.click();
    setTimeout(() => {
      URL.revokeObjectURL(url);
      a.remove();
    }, 0);
  };

  // Render
  return (
    <div className="w-full min-h-screen flex flex-col items-center gap-6 p-6">
      <h1 className="text-2xl font-semibold">
        3×3 Metaballs Grid — SVG Membrane
      </h1>

      <div className="grid grid-cols-6 gap-6 w-full max-w-5xl">
        <label className="flex flex-col gap-2">
          <span className="text-sm font-medium">Radius: {radius}px</span>
          <input
            type="range"
            min={8}
            max={80}
            step={1}
            value={radius}
            onChange={(e) => setRadius(parseInt(e.target.value, 10))}
            className="w-full accent-black"
          />
        </label>
        <label className="flex flex-col gap-2">
          <span className="text-sm font-medium">
            Connection reach factor: {reach.toFixed(2)}×
          </span>
          <input
            type="range"
            min={0.8}
            max={3}
            step={0.05}
            value={reach}
            onChange={(e) => setReach(parseFloat(e.target.value))}
            className="w-full accent-black"
          />
        </label>
        <label className="flex flex-col gap-2">
          <span className="text-sm font-medium">
            Gap between circles: {gap}px
          </span>
          <input
            type="range"
            min={0}
            max={120}
            step={1}
            value={gap}
            onChange={(e) => setGap(parseInt(e.target.value, 10))}
            className="w-full accent-black"
          />
        </label>
        <label className="flex flex-col gap-2">
          <span className="text-sm font-medium">
            Min metaball thickness (spread v): {vMin.toFixed(2)}
          </span>
          <input
            type="range"
            min={0}
            max={1}
            step={0.01}
            value={vMin}
            onChange={(e) => setVMin(parseFloat(e.target.value))}
            className="w-full accent-black"
          />
        </label>
        <label className="flex flex-col gap-2">
          <span className="text-sm font-medium">
            Stroke thickness: {strokeW}px
          </span>
          <input
            type="range"
            min={1}
            max={16}
            step={1}
            value={strokeW}
            onChange={(e) => setStrokeW(parseInt(e.target.value, 10))}
            className="w-full accent-black"
          />
        </label>
        <label className="flex flex-col gap-2">
          <span className="text-sm font-medium">Path color</span>
          <input
            type="color"
            value={pathColor}
            onChange={(e) => setPathColor(e.target.value)}
            className="h-10 w-full p-0 border rounded"
          />
        </label>
        <label className="flex flex-col gap-2">
          <span className="text-sm font-medium">Background color</span>
          <input
            type="color"
            value={bgColor}
            onChange={(e) => setBgColor(e.target.value)}
            className="h-10 w-full p-0 border rounded"
          />
        </label>
        <label className="flex flex-col gap-2">
          <span className="text-sm font-medium">
            Connection reach factor: {reach.toFixed(2)}× {autoReach ? "(auto)" : ""}
          </span>
          <input
            type="range"
            min={0.8}
            max={3}
            step={0.05}
            value={reach}
            onChange={(e) => setReach(parseFloat(e.target.value))}
            className="w-full accent-black disabled:opacity-50"
            disabled={autoReach}
          />
        </label>
        <label className="flex items-center gap-3">
          <input
            type="checkbox"
            checked={autoReach}
            onChange={(e) => setAutoReach(e.target.checked)}
            className="accent-black"
          />
          <span className="text-sm font-medium">Auto reach (recommended)</span>
        </label>
      </div>

      <div className="w-full max-w-5xl flex items-center justify-between">
        <button
          className="text-sm underline underline-offset-4"
          onClick={() => setAdvOpen((v) => !v)}
        >
          {advOpen ? "Hide Advanced" : "Show Advanced"}
        </button>
        <div className="flex items-center gap-3">
          <label className="flex items-center gap-2 text-sm">
            <input
              type="checkbox"
              checked={advAuto}
              onChange={(e) => setAdvAuto(e.target.checked)}
              className="accent-black"
            />
            Auto tolerances
          </label>
        </div>
      </div>

      {advOpen && (
        <div className="grid grid-cols-6 gap-6 w-full max-w-5xl">
          <label className="flex flex-col gap-2">
            <span className="text-sm font-medium">Edge probe (delta): {advDelta.toFixed(2)}px</span>
            <input
              type="range"
              min={0.2}
              max={10}
              step={0.1}
              value={advDelta}
              onChange={(e) => setAdvDelta(parseFloat(e.target.value))}
              disabled={advAuto}
              className="w-full accent-black disabled:opacity-50"
            />
          </label>
          <label className="flex flex-col gap-2">
            <span className="text-sm font-medium">Fill overlap: {advOverlap.toFixed(2)}px</span>
            <input
              type="range"
              min={0}
              max={12}
              step={0.1}
              value={advOverlap}
              onChange={(e) => setAdvOverlap(parseFloat(e.target.value))}
              disabled={advAuto}
              className="w-full accent-black disabled:opacity-50"
            />
          </label>
          <label className="flex flex-col gap-2">
            <span className="text-sm font-medium">Stitch gap: {advJoin.toFixed(2)}px</span>
            <input
              type="range"
              min={0}
              max={30}
              step={0.1}
              value={advJoin}
              onChange={(e) => setAdvJoin(parseFloat(e.target.value))}
              disabled={advAuto}
              className="w-full accent-black disabled:opacity-50"
            />
          </label>
          <label className="flex flex-col gap-2">
            <span className="text-sm font-medium">Min fragment: {advMinFrag.toFixed(2)}px</span>
            <input
              type="range"
              min={0}
              max={60}
              step={0.5}
              value={advMinFrag}
              onChange={(e) => setAdvMinFrag(parseFloat(e.target.value))}
              disabled={advAuto}
              className="w-full accent-black disabled:opacity-50"
            />
          </label>
        </div>
      )}

      <div className="w-full flex items-center justify-between max-w-5xl">
        <div />
        <button
          onClick={downloadSVG}
          className="px-3 py-2 rounded-xl border border-white/20 shadow-sm hover:shadow transition text-sm"
        >
          Download SVG
        </button>
      </div>

      <div className="w-full flex justify-center">
        <svg
          ref={svgRef}
          viewBox={`0 0 ${width} ${height}`}
          width={Math.min(900, width)}
          height={Math.min(900, height)}
          className="shadow rounded-2xl"
        >
          {/* Background */}
          <rect x={0} y={0} width={width} height={height} fill={bgColor} />

          {/* Pair blob hit-areas + fills (clickable inside shape at all times) */}
          {pairArtifacts.map((p, i) => (
            <path
              key={`fillpair-${i}`}
              d={p.fillD}
              fill={pairStates[i] === 2 ? pathColor : "transparent"}
              stroke="none"
              onClick={() => togglePair(i)}
              style={{ cursor: "pointer", pointerEvents: "all" }}
            />
          ))}

          {/* Circle hit-areas + fills (clickable inside circle at all times) */}
          {Array.from({ length: 3 }).map((_, r) =>
            Array.from({ length: 3 }).map((_, c) => {
              const key = `${r + 1},${c + 1}`;
              const [x, y] = centers[key];
              return (
                <circle
                  key={`fillcircle-${key}`}
                  cx={x}
                  cy={y}
                  r={radius}
                  fill={circleStates[key] === 2 ? pathColor : "transparent"}
                  stroke="none"
                  onClick={() => toggleCircle(key)}
                  style={{ cursor: "pointer", pointerEvents: "all" }}
                />
              );
            })
          )}

          {/* Membrane connectors (trimmed segments) */}
          {pairArtifacts.map((p, i) => (
            <g
              key={`pair-${i}`}
              onClick={() => togglePair(i)}
              style={{ cursor: "pointer" }}
            >
              {p.segments.map((d, j) => (
                <path
                  key={`seg-${i}-${j}`}
                  d={d}
                  fill="none"
                  stroke={pairStates[i] === 1 ? bgColor : pathColor}
                  strokeWidth={strokeW}
                  strokeLinecap="round"
                  strokeLinejoin="round"
                />
              ))}
            </g>
          ))}

          {/* Circle arcs / outlines */}
          {Array.from({ length: 3 }).map((_, r) =>
            Array.from({ length: 3 }).map((_, c) => {
              const key = `${r + 1},${c + 1}`;
              const arcs = arcsByCircle[key];
              const [x, y] = centers[key];
              const stroke = circleStates[key] === 1 ? bgColor : pathColor;
              return (
                <g
                  key={`circle-${key}`}
                  style={{ cursor: "pointer" }}
                  onClick={() => toggleCircle(key)}
                >
                  {arcs ? (
                    arcs.map((d, i) => (
                      <path
                        key={`arc-${key}-${i}`}
                        d={d}
                        fill="none"
                        stroke={stroke}
                        strokeWidth={strokeW}
                      strokeLinecap="round"
                      strokeLinejoin="round"
                      />
                  ))
                ) : (
                  <circle
                    cx={x}
                    cy={y}
                    r={radius}
                    fill="none"
                    stroke={stroke}
                    strokeWidth={strokeW}
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                )}
              </g>
            );
          })
          )}
        </svg>
      </div>

      <p className="text-xs text-neutral-300 max-w-3xl text-center">
        The connectors use the membrane method (no blurs/filters). Circles and
        membranes can be toggled between outline, hidden-outline, and filled
        states (click any circle edge or membrane). Colors and stroke width are
        exported with the "Download SVG" button.
      </p>

      <div className="text-xs text-neutral-400">
        Groups: G1 → (1,1) &amp; (2,2). G2 → (1,2), (2,3), (3,2). Only
        intra-group connections are drawn.
      </div>
    </div>
  );
}

// ---------------- Dev Unit Tests (lightweight, run once) ----------------
// NOTE: These are console-based sanity checks — no UI changes and no dependencies.
(function runUnitTestsOnce() {
  if (typeof window === "undefined") return;
  const w = window as any;
  if (w.__MB_TESTS_RAN__) return; // avoid re-running on HMR
  w.__MB_TESTS_RAN__ = true;

  const approx = (a: number, b: number, tol = 1e-6) => Math.abs(a - b) <= tol;
  const log = (...args: any[]) => console.log("[MetaballsGrid tests]", ...args);

  try {
    // Test 1: XML header is properly terminated with newline
    console.assert(
      XML_HEADER.endsWith("\n"),
      "XML header should end with newline"
    );
    console.assert(
      XML_HEADER.includes("<?xml"),
      "XML header should include xml declaration"
    );

    // Test 2: minorSpan never exceeds PI
    const [s1, e1] = minorSpan(0.1, 5.9);
    console.assert(
      e1 - s1 <= Math.PI + 1e-9,
      "minorSpan length should be <= PI"
    );

    // Test 3: complementSegments basic overlap case
    const comp = complementSegments(
      mergeSpans([
        [0.1, 1.1],
        [1.0, 2.0],
      ])
    );
    // Expected: two segments: [0, 0.1] and [2.0, TAU]
    console.assert(
      comp.length === 2,
      "complementSegments should produce two ranges"
    );
    console.assert(
      approx(comp[0][0], 0) && approx(comp[0][1], 0.1),
      "first range should be [0, 0.1]"
    );
    console.assert(
      approx(comp[1][0], 2.0) && approx(comp[1][1], TAU),
      "second range should be [2.0, TAU]"
    );

    // Test 4: metaballParams returns null for far circles and non-null for near ones
    const far = metaballParams(30, 30, [0, 0], [1000, 0], 2.4, 0.5);
    console.assert(
      far === null,
      "metaballParams should be null for very distant circles"
    );
    const near = metaballParams(30, 30, [0, 0], [50, 0], 2.4, 0.5);
    console.assert(
      near !== null,
      "metaballParams should return params for reasonable distances"
    );

    log("All tests passed.");
  } catch (err) {
    console.warn("Unit tests encountered an error:", err);
  }
})();
