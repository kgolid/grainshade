import SimplexNoise from 'simplex-noise';

const simplex = new SimplexNoise();

export function sum_octave(num_iterations, x, y, scale, persistence, sigmoid_intensity) {
  let noise = 0;
  let maxAmp = 0;
  let amp = 1;
  let freq = 1 / scale;

  for (let i = 0; i < num_iterations; i++) {
    noise += simplex.noise3D(x * freq, y * freq, i) * amp;
    maxAmp += amp;
    amp *= persistence;
    freq *= 2;
  }
  var output = apply_sigmoid(noise / maxAmp, sigmoid_intensity);
  return output;
}

function apply_sigmoid(value, intensity) {
  if (intensity === 0) return value;
  return 2 * sigmoid(value * intensity) - 1;
}

function sigmoid(x) {
  return 1 / (1 + Math.exp(-x));
}
