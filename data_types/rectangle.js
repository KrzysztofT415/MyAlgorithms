export class Rectangle {
    constructor(x, y, w, h) {
        this.x = x
        this.y = y
        this.w = w
        this.h = h
    }

    *[Symbol.iterator]() {
        yield this.x
        yield this.y
        yield this.w
        yield this.h
    }

    contains(point) {
        return (
            this.x <= point.x && //
            point.x <= this.x + this.w &&
            this.y <= point.y &&
            point.y <= this.y + this.h
        )
    }

    intersects(other) {
        return !(
            other.x > this.x + this.w || //
            other.x + other.w < this.x ||
            other.y > this.y + this.h ||
            other.y + other.h < this.y
        )
    }
}
