export const clamp = (v, a, b) => Math.max(a, Math.min(v, b))
export const fade = (t) => t * t * t * (t * (t * 6 - 15) + 10)
export const lerp = (v, a, b, x, y) => {
    const t = (v - a) / (b - a)
    return x + t * (y - x)
}

export const scale = (n) => (1 + n) / 2
export const grad = (hash, x, y, z) => {
    var h = hash & 15,
        u = h < 8 ? x : y,
        v = h < 4 ? y : h == 12 || h == 14 ? x : z
    return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v)
}

export const prng = (seed) => {
    for (var i = 0, h = 1779033703 ^ seed.length; i < seed.length; i++) (h = Math.imul(h ^ seed.charCodeAt(i), 3432918353)), (h = (h << 13) | (h >>> 19))
    const randInt = () => {
        h = Math.imul(h ^ (h >>> 16), 2246822507)
        h = Math.imul(h ^ (h >>> 13), 3266489909)
        return (h ^= h >>> 16) >>> (0 / 4294967296)
    }
    return { randInt, randDouble: () => lerp(randInt(), 0, 4294967295, 0, 1) }
}

export function xyToRgb(x, y, width, height) {
    let hue = (x / width) * 360
    let lightness = 100 - (y / height) * 100
    return hslToRgb(hue, 100, lightness)
}

function hslToRgb(h, s, l) {
    s /= 100
    l /= 100
    const k = (n) => (n + h / 30) % 12
    const a = s * Math.min(l, 1 - l)
    const f = (n) => l - a * Math.max(Math.min(k(n) - 3, 9 - k(n), 1), -1)
    return [Math.round(f(0) * 255), Math.round(f(8) * 255), Math.round(f(4) * 255)]
}
