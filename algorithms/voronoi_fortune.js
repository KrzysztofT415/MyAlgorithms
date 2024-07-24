import { PriorityQueue } from '../data_structures/priorityqueue.js'

let opposites = { l: 'r', r: 'l' }

class Geometry {
    // Returns an edge object that represents the line (in the form of ax + by = c) that bisects two sites.
    static bisect(s1, s2) {
        let m = (s2.x - s1.x) / (s1.y - s2.y)
        let [mx, my] = [Math.abs((s1.x + s2.x) / 2), Math.abs((s1.y + s2.y) / 2)]
        let [a, b, c] = [-1 * m, 1, -1 * (m * mx - my)]
        return { region: { l: s1, r: s2 }, ep: { l: null, r: null }, a, b, c } // prettier-ignore
    }

    // Returns a coordinate (x,y) that is the intersection of two half edge objects.
    static intersect(el1, el2) {
        let [e1, e2] = [el1.edge, el2.edge]
        if (!e1 || !e2 || e1.region.r == e2.region.r) return null

        let d = e1.a * e2.b - e1.b * e2.a
        if (Math.abs(d) < 1e-10) return null

        let [xint, yint, e1r, e2r] = [(e1.c * e2.b - e2.c * e1.b) / d, (e2.c * e1.a - e1.c * e2.a) / d, e1.region.r, e2.region.r]

        let [el, e] = [el2, e2]
        if (e1r.y < e2r.y || (e1r.y == e2r.y && e1r.x < e2r.x)) [el, e] = [el1, e1]

        let rightOfSite = xint >= e.region.r.x

        if ((rightOfSite && el.side === 'l') || (!rightOfSite && el.side === 'r')) return null

        return {
            x: xint,
            y: yint,
        }
    }

    // Returns a bool.
    static rightOf(he, p) {
        let e = he.edge,
            topsite = e.region.r,
            rightOfSite = p.x > topsite.x

        if (rightOfSite && he.side === 'l') return 1
        if (!rightOfSite && he.side === 'r') return 0

        let yl = e.c - e.a * p.x,
            t1 = p.y - yl,
            t2 = p.x - topsite.x,
            t3 = yl - topsite.y

        let above = t1 * t1 > t2 * t2 + t3 * t3
        return he.side === 'l' ? above : !above
    }

    // Used to reconstruct the diagram.
    static endPoint(edge, side, site, callback) {
        edge.ep[side] = site
        if (!edge.ep[opposites[side]]) return
        callback(edge)
    }

    // Computes the Euclidean distance between two coordinates (x, y)
    static distance(s, t, type) {
        let dx = s.x - t.x,
            dy = s.y - t.y
        if (type == 'm') return Math.abs(dx) + Math.abs(dy) // manhattan
        if (type == 's') return Math.pow(dx * dx * dx + dy * dy * dy, 1 / 3) // smooth
        return Math.sqrt(dx * dx + dy * dy) // euclidean
    }

    static clipPolygonToBoundary(subjectPolygon, boundary) {
        let outputList = subjectPolygon
        let cp1, cp2, s, e

        const clipPolygon = [
            [boundary.x, boundary.y],
            [boundary.x + boundary.w, boundary.y],
            [boundary.x + boundary.w, boundary.y + boundary.h],
            [boundary.x, boundary.y + boundary.h],
        ]

        for (let j = 0; j < clipPolygon.length; j++) {
            const inputList = outputList
            outputList = []
            cp1 = clipPolygon[j]
            cp2 = clipPolygon[(j + 1) % clipPolygon.length]

            for (let i = 0; i < inputList.length; i++) {
                s = inputList[i]
                e = inputList[(i + 1) % inputList.length]

                if (Geometry.inside(e, cp1, cp2)) {
                    if (!Geometry.inside(s, cp1, cp2)) {
                        outputList.push(Geometry.intersection(cp1, cp2, s, e))
                    }
                    outputList.push(e)
                } else if (Geometry.inside(s, cp1, cp2)) {
                    outputList.push(Geometry.intersection(cp1, cp2, s, e))
                }
            }
        }
        return outputList
    }

    static inside(p, cp1, cp2) {
        return (cp2[0] - cp1[0]) * (p[1] - cp1[1]) > (cp2[1] - cp1[1]) * (p[0] - cp1[0])
    }

    static intersection(cp1, cp2, s, e) {
        const dc = [cp1[0] - cp2[0], cp1[1] - cp2[1]]
        const dp = [s[0] - e[0], s[1] - e[1]]
        const n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0]
        const n2 = s[0] * e[1] - s[1] * e[0]
        const n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0])
        return [(n1 * dp[0] - n2 * dc[0]) * n3, (n1 * dp[1] - n2 * dc[1]) * n3]
    }
}

class EdgeList {
    constructor() {
        this.leftEnd = this.createHalfEdge(null, 'l')
        this.rightEnd = this.createHalfEdge(null, 'l')
        this.leftEnd.r = this.rightEnd
        this.rightEnd.l = this.leftEnd
    }

    // Constructs a half edge object.
    createHalfEdge(edge, side) {
        return {
            edge: edge,
            side: side,
            vertex: null,
            l: null,
            r: null,
        }
    }

    // he represents a half edge.
    insert(lb, he) {
        he.l = lb
        he.r = lb.r
        lb.r.l = he
        lb.r = he
    }

    // site: site object
    // returns a half edge (leaf node?)
    leftBound(site) {
        let he = this.leftEnd
        he = he.r
        while (he != this.rightEnd && Geometry.rightOf(he, site)) he = he.r
        he = he.l
        return he
    }

    del(he) {
        he.l.r = he.r
        he.r.l = he.l
        he.edge = null
    }

    leftRegion(he) {
        if (he.edge == null) return null
        return he.edge.region[he.side]
    }

    rightRegion(he) {
        if (he.edge == null) return null
        return he.edge.region[opposites[he.side]]
    }
}

export default class VoronoiDiagram {
    constructor(width, height) {
        this.width = width
        this.height = height
        this.sites = []
        this.voronoi = []
        this.rendersImage = false
    }

    pointToPriority(p) {
        return p.ystar * this.width ** 2 + p.vertex.x
    }

    twoPointsCase() {
        let [a, b, c, d] = [this.sites[0].x, this.sites[0].y, this.sites[1].x, this.sites[1].y]
        let [width, height] = [this.width, this.height]
        const midX = (a + c) / 2
        const midY = (b + d) / 2
        const slope = (d - b) / (c - a)
        const perpendicularSlope = -1 / slope

        const lineIntersection = (x1, y1, x2, y2) => {
            const intercept = midY - perpendicularSlope * midX
            let points = []

            let y = intercept
            if (y >= 0 && y <= height) points.push([0, y])
            y = perpendicularSlope * width + intercept
            if (y >= 0 && y <= height) points.push([width, y])
            let x = (0 - intercept) / perpendicularSlope
            if (x >= 0 && x <= width) points.push([x, 0])
            x = (height - intercept) / perpendicularSlope
            if (x >= 0 && x <= width) points.push([x, height])

            return points
        }

        const [p1, p2] = lineIntersection(0, 0, width, height)
        const isHorizontal = Math.abs(slope) < 1
        let polygon1, polygon2

        if (!isHorizontal) {
            if (p2.x == 0) [p1, p2] = [p2, p1]
            polygon1 = { borders: [[0, 0], p1, p2, [width, 0]], center: this.sites[0] }
            polygon2 = { borders: [p1, [0, height], [width, height], p2], center: this.sites[1] }
            if (d > b) this.voronoi = [polygon1, polygon2]
            else this.voronoi = [polygon2, polygon1]
        } else {
            if (p2.y == 0) [p1, p2] = [p2, p1]
            polygon1 = { borders: [[0, 0], p1, p2, [0, height]], center: this.sites[0] }
            polygon2 = { borders: [p1, [width, 0], [width, height], p2], center: this.sites[1] }
            if (c > a) this.voronoi = [polygon1, polygon2]
            else this.voronoi = [polygon2, polygon1]
        }
    }

    calculateVoronoi() {
        if (this.sites.length == 1) {
            this.voronoi = [{borders: [[0, 0], [this.width, 0], [this.width, this.height], [0, this.height],], center: this.sites[0]}] // prettier-ignore
            return
        }

        if (this.sites.length == 2) {
            // TODO: fix proper dividing when 2 points
            this.twoPointsCase()
            return
        }

        let polygons = this.sites.map(() => [])

        const callback = (e) => {
            let s1, s2, x1, x2, y1, y2

            if (e.a === 1 && e.b >= 0) {
                s1 = e.ep.r
                s2 = e.ep.l
            } else {
                s1 = e.ep.l
                s2 = e.ep.r
            }
            if (e.a === 1) {
                y1 = s1 ? s1.y : -1e6
                x1 = e.c - e.b * y1
                y2 = s2 ? s2.y : 1e6
                x2 = e.c - e.b * y2
            } else {
                x1 = s1 ? s1.x : -1e6
                y1 = e.c - e.a * x1
                x2 = s2 ? s2.x : 1e6
                y2 = e.c - e.a * x2
            }
            let v1 = [x1, y1],
                v2 = [x2, y2]

            polygons[e.region.l.index].push(v1, v2)
            polygons[e.region.r.index].push(v1, v2)
        }

        let sitesList = this.sites
            .map((v, i) => {
                return {
                    index: i,
                    x: v.x,
                    y: v.y,
                }
            })
            .sort((a, b) => (a.y < b.y ? -1 : a.y > b.y ? 1 : a.x < b.x ? -1 : a.x > b.x ? 1 : 0))

        const edgeList = new EdgeList()
        const eventQueue = new PriorityQueue()

        let firstSite = sitesList.shift()
        let newSite = sitesList.shift()
        let newIntStar
        let lbnd, rbnd, llbnd, rrbnd, bisector
        let bot, top, p, v
        let e, pm

        while (true) {
            if (!eventQueue.is_empty()) {
                let elem = eventQueue.front
                newIntStar = { x: elem.vertex.x, y: elem.ystar }
            }

            if (newSite && (eventQueue.is_empty() || newSite.y < newIntStar.y || (newSite.y == newIntStar.y && newSite.x < newIntStar.x))) {
                lbnd = edgeList.leftBound(newSite)
                rbnd = lbnd.r

                bot = edgeList.rightRegion(lbnd) ?? firstSite
                e = Geometry.bisect(bot, newSite)

                bisector = edgeList.createHalfEdge(e, 'l')
                edgeList.insert(lbnd, bisector)
                p = Geometry.intersect(lbnd, bisector)
                if (p) {
                    lbnd.vertex = p
                    lbnd.ystar = p.y + Geometry.distance(p, newSite)
                    eventQueue.update(lbnd, this.pointToPriority(lbnd))
                }

                lbnd = bisector
                bisector = edgeList.createHalfEdge(e, 'r')
                edgeList.insert(lbnd, bisector)
                p = Geometry.intersect(bisector, rbnd)
                if (p) {
                    bisector.vertex = p
                    bisector.ystar = p.y + Geometry.distance(p, newSite)
                    eventQueue.insert(bisector, this.pointToPriority(bisector))
                }
                newSite = sitesList.shift()
            } else if (!eventQueue.is_empty()) {
                lbnd = eventQueue.pop()
                llbnd = lbnd.l
                rbnd = lbnd.r
                rrbnd = rbnd.r
                bot = edgeList.leftRegion(lbnd) ?? firstSite
                top = edgeList.rightRegion(rbnd, firstSite) ?? firstSite
                v = lbnd.vertex

                Geometry.endPoint(lbnd.edge, lbnd.side, v, callback)
                Geometry.endPoint(rbnd.edge, rbnd.side, v, callback)
                edgeList.del(lbnd)
                eventQueue.remove(rbnd)
                edgeList.del(rbnd)

                pm = 'l'
                if (bot.y > top.y) [bot, top, pm] = [top, bot, 'r']

                e = Geometry.bisect(bot, top)
                bisector = edgeList.createHalfEdge(e, pm)

                edgeList.insert(llbnd, bisector)
                Geometry.endPoint(e, opposites[pm], v, callback)

                p = Geometry.intersect(llbnd, bisector)
                if (p) {
                    llbnd.vertex = p
                    llbnd.ystar = p.y + Geometry.distance(p, bot)
                    eventQueue.update(llbnd, this.pointToPriority(llbnd))
                }
                p = Geometry.intersect(bisector, rrbnd)
                if (p) {
                    bisector.vertex = p
                    bisector.ystar = p.y + Geometry.distance(p, bot)
                    eventQueue.insert(bisector, this.pointToPriority(bisector))
                }
            } else break
        }

        for (lbnd = edgeList.leftEnd.r; lbnd != edgeList.rightEnd; lbnd = lbnd.r) callback(lbnd.edge)

        // Reconnect the polygon segments into counterclockwise loops.
        polygons = polygons.map((polygon, index) => {
            let [cx, cy] = [this.sites[index].x, this.sites[index].y]
            let [w2, h2] = [this.width ** 3, this.height ** 3]
            polygon.forEach((v) => {
                v.angle = Math.atan2(v[0] - cx, v[1] - cy)
                v[0] = Math.min(Math.max(-w2, v[0]), w2)
                v[1] = Math.min(Math.max(-h2, v[1]), h2)
            })
            polygon = polygon
                .sort((a, b) => a.angle - b.angle) //
                .filter((d, i) => !i || d.angle - polygon[i - 1].angle > 1e-10)

            return {
                borders: Geometry.clipPolygonToBoundary(polygon, this.sites[index].boundary),
                center: this.sites[index],
            }
        })

        this.voronoi = polygons
    }
}
