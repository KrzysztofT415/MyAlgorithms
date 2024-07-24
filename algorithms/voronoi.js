export default class VoronoiDiagram {
    constructor(width, height) {
        this.width = width
        this.height = height
        this.sites = []
        this.voronoi = []
        this.rendersImage = true
        this.distance_type = 'manhattan'
        this.degree = 3
    }

    distance(x1, x2, y1, y2) {
        if (this.distance_type == 'euclidean') return Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) // Euclidean distance
        if (this.distance_type == 'manhattan') return Math.abs(x1 - x2) + Math.abs(y1 - y2) // Manhattan distance
        return Math.pow(Math.pow(Math.abs(x1 - x2), this.degree) + Math.pow(Math.abs(y1 - y2), this.degree), 1 / this.degree) // Smooth distance
    }

    sample(x, y) {
        const neighbours = this.getNeighbours(x, y)
        if (neighbours.length == 0) return 0
        let n = 0
        let d = this.distance(neighbours[0].x, x, neighbours[0].y, y)
        neighbours.forEach((neighbour, index) => {
            const d_tmp = this.distance(neighbour.x, x, neighbour.y, y)
            if (d_tmp < d) [n, d] = [index, d_tmp]
        })
        return this.sites.indexOf(neighbours[n])
    }

    makeMap() {
        const map = Array.from(Array(this.width), () => new Array(this.height))
        for (let y = 0; y < this.height; y++) {
            for (let x = 0; x < this.width; x++) {
                map[x][y] = this.sample(x, y)
            }
        }
        return map
    }

    makeImageData() {
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

    calculateVoronoi(getNeighbours) {
        this.getNeighbours = getNeighbours ?? (() => this.sites)
        if (this.sites.length == 0) return
        this.voronoi = this.makeImageData()
    }
}
