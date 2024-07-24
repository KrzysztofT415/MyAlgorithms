export class Particle {
    constructor(x, y, v, r, boundary) {
        this.x = x
        this.y = y
        this.v = v
        this.r = r
        this.boundary = boundary
    }

    *[Symbol.iterator]() {
        yield this.x
        yield this.y
        yield this.v
        yield this.boundary
    }

    update() {
        if (!this.v.x && !this.v.y) return
        if (this.x + this.v.x < this.boundary.x || this.boundary.w < this.x + this.v.x) this.v.x *= -1
        if (this.y + this.v.y < this.boundary.y || this.boundary.h < this.y + this.v.y) this.v.y *= -1
        this.x += this.v.x
        this.y += this.v.y
    }
}
