import PoissonDiskSampling from 'poisson-disk-sampling';

import { sum_octave } from './noise';

let sketch = function (p) {
  let THE_SEED;
  let noise_grid;
  let slope_grid;
  let pds;
  let points;

  const darkness = 8;

  const scale = 500;
  const persistence = 0.35;
  const sigmoid_intensity = 4;

  const grid_dim_x = 1000;
  const grid_dim_y = 1000;
  const padding = 40;
  const canvas_dim_x = grid_dim_x + 2 * padding;
  const canvas_dim_y = grid_dim_y + 2 * padding;
  const cell_dim = 2;
  const nx = grid_dim_x / cell_dim;
  const ny = grid_dim_y / cell_dim;

  p.setup = function () {
    p.createCanvas(canvas_dim_x, canvas_dim_y);
    p.pixelDensity(4);
    p.noLoop();

    THE_SEED = p.floor(p.random(9999999));
    p.randomSeed(THE_SEED);
  };

  p.draw = function () {
    p.background('#fd0');
    reset();
    displaySample();
  };

  p.keyPressed = function () {
    if (p.keyCode === 80) p.saveCanvas('sketch_' + THE_SEED, 'jpeg');
  };

  function displaySample() {
    p.stroke(0);
    p.translate(padding, padding);
    for (const pnt of points) {
      p.strokeWeight(0.1 + Math.random() * 2.5);
      p.point(pnt[0], pnt[1]);
    }
  }

  function reset() {
    noise_grid = build_noise_grid(0.5, 0);
    slope_grid = build_slope_grid();

    pds = createPoisson(slope_grid);
    points = pds.fill();

    p.randomSeed(THE_SEED);
  }

  function build_noise_grid(baseline, offset_mag) {
    return [...Array(ny + 1)].map((_, y) =>
      [...Array(nx + 1)].map(
        (_, x) =>
          sum_octave(16, x, y, scale, persistence, sigmoid_intensity) +
          offset_mag * (center_offset(x, y) - baseline)
      )
    );
  }

  function build_slope_grid() {
    const grid = [];
    for (let i = 0; i < noise_grid.length; i++) {
      const row = [];
      for (let j = 0; j < noise_grid[i].length; j++) {
        row.push(map_to_diff(j, i));
      }
      grid.push(row);
    }
    return grid;
  }

  function createPoisson(slope) {
    return new PoissonDiskSampling({
      shape: [grid_dim_x, grid_dim_y],
      minDistance: 1,
      maxDistance: 15,
      distanceFunction: function (p) {
        return Math.pow(slope[Math.floor(p[1] / cell_dim)][Math.floor(p[0] / cell_dim)], 8); // value between 0 and 1
      },
    });
  }

  // Output range [-1, 1]
  function center_offset(x, y) {
    return 1 - distance_from_centre(x, y) * 3;
  }

  // Output range: [0, 1] (within tangent circle);
  function distance_from_centre(x, y) {
    return Math.sqrt(Math.pow(nx / 2 - x, 2) + Math.pow(ny / 2 - y, 2)) / nx;
  }

  function map_to_diff(x, y) {
    const slope = get_slope(x, y);
    const mag = Math.sqrt(Math.pow(slope[0], 2) + Math.pow(slope[1], 2));
    const dir = get_dir(...slope);
    const dirdiff = Math.abs(2 + dir);
    return p.constrain(1 - mag * dirdiff * darkness, 0, 1);
  }

  function get_slope(x, y) {
    if (y <= 0 || y >= noise_grid.length - 1) return [0, 0];
    if (x <= 0 || x >= noise_grid[y].length - 1) return [0, 0];
    const n = noise_grid[y - 1][x];
    const s = noise_grid[y + 1][x];
    const w = noise_grid[y][x - 1];
    const e = noise_grid[y][x + 1];

    const lat = s - n;
    const lon = e - w;

    return [lon, lat];
  }

  function get_dir(w, h) {
    let v = p.createVector(w, h);
    return v.heading();
  }
};
new p5(sketch);
