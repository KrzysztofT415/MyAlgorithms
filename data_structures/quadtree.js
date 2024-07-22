import { Rectangle } from '../data_types/rectangle.js'

class QuadTreeData {
    constructor(x, y, data) {
        this.x = x
        this.y = y
        this.data_stored = data
        this.owner = null
    }
}

class QuadTreeNode {
    constructor(parent, boundary, c) {
        this.parent = parent
        this.boundary = boundary
        this.capacity = c
        this.data_stored = []
        this.children = []
    }

    *getData() {
        if (this.isLeaf()) for (const data of this.data_stored) yield data
        for (const child of this.children) yield* child.getData()
    }

    *getNodes() {
        if (this.isLeaf()) yield this
        for (const child of this.children) yield* child.getNodes()
    }

    contains = (p) => this.boundary.contains(p)
    isLeaf = () => this.children.length == 0

    insert = (data) => {
        if (this.isLeaf()) {
            data.owner = this
            this.data_stored.push(data)
            if (this.data_stored.length > this.capacity) this.subdivide()
            return
        }
        for (const child of this.children) {
            if (child.contains(data)) {
                child.insert(data)
                return
            }
        }
    }

    subdivide = () => {
        const [x, y, w, h, c] = [this.boundary.x, this.boundary.y, this.boundary.w / 2, this.boundary.h / 2, this.capacity]
        this.children[0] = new QuadTreeNode(this, new Rectangle(x, y, w, h), c)
        this.children[1] = new QuadTreeNode(this, new Rectangle(x + w, y, w, h), c)
        this.children[2] = new QuadTreeNode(this, new Rectangle(x, y + h, w, h), c)
        this.children[3] = new QuadTreeNode(this, new Rectangle(x + w, y + h, w, h), c)

        for (const data of this.data_stored) this.insert(data)
        this.data_stored = []
    }

    remove = (data) => {
        if (this.isLeaf()) {
            this.data_stored = this.data_stored.filter((p) => p.x != data.x && p.y != data.y)
            if (this.parent.getDataLength() <= this.capacity) this.parent.reverseSubdivision()
            return
        }
        for (const child of this.children) {
            if (child.contains(data)) {
                child.remove(data)
                return
            }
        }
    }

    getDataLength = () => {
        if (this.isLeaf()) return this.data_stored.length
        let sum = 0
        for (const child of this.children) sum += child.getDataLength()
        return sum
    }

    reverseSubdivision = () => {
        if (this.isLeaf()) {
            this.parent.children = []
            for (const data of this.data_stored) this.parent.insert(data)
            this.data_stored = []
        }
        for (const child of this.children) child.reverseSubdivision()
    }

    getNode = (data) => {
        if (this.isLeaf()) return this
        for (const child of this.children) {
            if (child.contains(data)) {
                return child.getNode(data)
            }
        }
    }

    query = (range, found) => {
        if (!range.intersects(this.boundary)) return found
        if (this.isLeaf()) {
            for (let data of this.data_stored) {
                if (range.contains(data)) found.push(data.data_stored)
            }
            return found
        }
        for (const child of this.children) {
            child.query(range, found)
        }
        return found
    }
}

export class QuadTree {
    constructor(x, y, w, h, capacity) {
        this.root = new QuadTreeNode(undefined, new Rectangle(x, y, w, h), capacity)
    }

    *getData() {
        yield* this.root.getData()
    }

    *getNodes() {
        yield* this.root.getNodes()
    }

    insert = (x, y, data) => this.root.insert(new QuadTreeData(x, y, data))
    query = (range) => this.root.query(range, [])
    find = (x, y) => this.root.getNode({ x, y })
}
