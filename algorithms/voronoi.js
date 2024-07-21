import { xyToRgb } from '../utils.js'

export class VoronoiDiagram {
    constructor(width, height) {
        this.width = width
        this.height = height
        this.sites = []
        this.voronoi = []
        this.isPlane = true
    }

    distance(x1, x2, y1, y2) {
        if (this.params.type == 'euclidean') return Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) // Euclidean distance
        if (this.params.type == 'manhattan') return Math.abs(x1 - x2) + Math.abs(y1 - y2) // Manhattan distance
        return Math.pow(Math.pow(Math.abs(x1 - x2), this.params.degree) + Math.pow(Math.abs(y1 - y2), this.params.degree), 1 / this.params.degree) // Smooth distance
    }

    sample(x, y) {
        let n = 0
        let d = this.distance(this.sites[0].x, x, this.sites[0].y, y)
        for (let c = 1; c < this.sites.length; c++) {
            const d_tmp = this.distance(this.sites[c].x, x, this.sites[c].y, y)
            if (d_tmp < d) [n, d] = [c, d_tmp]
        }
        return n
    }

    makeMap = () => {
        const map = Array.from(Array(this.width), () => new Array(this.height))
        for (let y = 0; y < this.height; y++) {
            for (let x = 0; x < this.width; x++) {
                map[x][y] = this.sample(x, y)
            }
        }
        return map
    }

    makeImageData = () => {
        const data = []
        for (let y = 0; y < this.height; y++) {
            for (let x = 0; x < this.width; x++) {
                const i = (y * this.width + x) * 4
                const [r, g, b] = this.colors[this.sample(x, y)]
                data[i] = r
                data[i + 1] = g
                data[i + 2] = b
                data[i + 3] = 255
            }
        }
        return data
    }

    calculateVoronoi = () => {
        this.colors = this.sites.map((p) => xyToRgb(p.x, p.y, this.width, this.height))
        if (this.sites.length == 0) return
        this.imageData = this.makeImageData()
    }
}
