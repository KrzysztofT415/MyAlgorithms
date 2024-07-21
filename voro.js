// import { prngd, xyToRgb } from './utils.js'
// import { Point } from './point.js'
// import { PriorityQueue } from './libs/priorityqueue.js'
// import { doVoronoi } from './voro.js'

// // Just a binary tree. leaves are arcs. inner nodes are breakpoints.
// class BeachLine {
//     constructor() {
//         this.root = null
//     }

//     empty = () => this.root === null

//     initialInsert = (p) => (this.root = this.createArc(p, null))

//     /*
//      * brokenArc: a leaf node (an arc) of the beach line
//      * newSite: the new site event
//      * edge: the perpendicular bisector of brokenArc.p and (newSite.x, newSite.y)
//      *
//      * Splits a single arc into three arcs. This means a leaf node
//      * is transformed into a subtree containing two inner nodes (two new breakpoints).
//      */
//     split = (brokenArc, newSite, edge) => {
//         let newP = { x: newSite.x, y: newSite.y }
//         let lbp = this.createBp(brokenArc.p, newP, edge, brokenArc.side)
//         let rbp = this.createBp(newP, brokenArc.p, edge, 'r')
//         lbp.r = rbp
//         rbp.parent = lbp
//         lbp.l = this.createArc(brokenArc.p, 'l')
//         lbp.l.parent = lbp
//         let newArc = this.createArc(newP, 'l')
//         rbp.l = newArc
//         rbp.l.parent = rbp
//         rbp.r = this.createArc(brokenArc.p, 'r')
//         rbp.r.parent = rbp

//         if (brokenArc.parent) {
//             brokenArc.parent[brokenArc.side] = lbp
//             lbp.parent = brokenArc.parent
//         } else this.root = lbp

//         return {
//             l: lbp.l,
//             m: newArc,
//             r: rbp.r,
//         }
//     }

//     flip = (side) => (side === 'l' ? 'r' : 'l')
//     // Removes the arc deleteArc from the tree.
//     // Probably the most difficult part to understand.
//     remove = (deleteArc, y) => {
//         let lArc = this.leftArc(deleteArc)
//         let rArc = this.rightArc(deleteArc)

//         let grandParent
//         if (deleteArc.side === 'l') grandParent = this.commonAncestor(deleteArc, lArc)
//         else grandParent = this.commonAncestor(deleteArc, rArc)

//         grandParent.regions.l = lArc.p
//         grandParent.regions.r = rArc.p
//         grandParent.edge = Geometry.bisect(lArc.p, rArc.p)

//         let removed = deleteArc.parent

//         // The node opposite of deleteArc is being moved, so its side may need to be adjusted.
//         removed[this.flip(deleteArc.side)].side = removed.side

//         // Removes the breakpoint from the tree.
//         // Conceptually this line does this:
//         // If we have A <-> B <-> C, then A is changed to point at C (B is the letiable removed).
//         removed.parent[removed.side] = removed[this.flip(deleteArc.side)]
//         // But C still has a link to B, so point C at A.
//         removed[this.flip(deleteArc.side)].parent = removed.parent
//         // Now we have A <-> C, but C may have switched sides in the tree, so fix it.
//         // It will just be the side of removed.
//         removed[this.flip(deleteArc.side)].side = removed.side

//         return {
//             l: lArc,
//             r: rArc,
//         }
//     }

//     /*
//      * Returns the arc (a leaf in the tree) that is above the given location.
//      * In other terms, imagine a vertical ray shooting up from the site...
//      * which arc does it hit?
//      *
//      * This performs a binary search on the breakpoints (the internal nodes)
//      * to reach the leaf.
//      */
//     search = (site) => this.bpSearch(site, this.root)

//     // Logic for search
//     bpSearch = (site, node) => {
//         if (!node.l && !node.r) return node

//         let breakpoint = Geometry.breakpoint(node.regions.l, node.regions.r, node.edge, site.y)

//         if (breakpoint > site.x) return this.bpSearch(site, node.l)
//         else if (breakpoint < site.x) return this.bpSearch(site, node.r)
//     }

//     // Returns the arc right of arc or null if arc is the right-most arc.
//     // Go one leaf right.
//     rightArc = (arc) => {
//         let node = arc
//         while (node.side === 'r' || node.side === null) {
//             if (!node.parent) return null
//             node = node.parent
//         }
//         node = node.parent.r
//         while (node.type != 'arc') node = node.l
//         return node
//     }

//     // Returns the arc left of arc or null if arc is the left-most arc.
//     // Go one leaf left.
//     leftArc = (arc) => {
//         let node = arc
//         while (node.side === 'l' || node.side === null) {
//             if (!node.parent) return null
//             else node = node.parent
//         }
//         node = node.parent.l
//         while (node.type != 'arc') node = node.r
//         return node
//     }

//     // Returns the left-most arc.
//     leftMostArc = () => {
//         let node = root
//         while (node.type === 'breakpoint') node = node.l
//         return node
//     }

//     depth = (node) => {
//         let d = 0
//         while (node.parent) {
//             node = node.parent
//             d++
//         }
//         return d
//     }

//     // Find the first common ancestor of two nodes.
//     // Used in delete to rearrange tree.
//     commonAncestor(n1, n2) {
//         let d1 = this.depth(n1)
//         let d2 = this.depth(n2)

//         if (d1 === d2) return this.walkUpTree(n1, n2)

//         // Adjust depths to the same level.
//         let lower, higher, lowerNode, higherNode
//         d1 < d2 ? ((lower = d2), (lowerNode = n2), (higher = d1), (higherNode = n1)) : ((lower = d1), (lowerNode = n1), (higher = d2), (higherNode = n2))

//         while (lower > higher) {
//             lowerNode = lowerNode.parent
//             lower -= 1
//         }

//         return this.walkUpTree(lowerNode, higherNode)
//     }

//     // Helper to commonAncestor.
//     walkUpTree(n1, n2) {
//         while (n1 != n2) {
//             n1 = n1.parent
//             n2 = n2.parent
//         }
//         return n1 // or n2.
//     }

//     // Create a leaf node (arc) in the tree.
//     createArc = (p, side) => {
//         return {
//             p: p,
//             parent: null,
//             side: side,
//             type: 'arc',
//         }
//     }

//     // Create an inner node (breakpoint) in the tree.
//     createBp = (l, r, edge, side) => {
//         return {
//             regions: {
//                 l: l,
//                 r: r,
//             },
//             edge: edge,
//             side: side,
//             parent: null,
//             type: 'breakpoint',
//         }
//     }
// }

// // Mathematical operations.
// class Geometry {
//     /*
//      * Returns the perpendicular bisector of two points in the form of ax + by = c as the triplet (a,b,c).
//      */
//     static bisect = (s1, s2) => {
//         let m = Geometry.reciprocal(Geometry.slope(s1, s2))
//         let midpoint = Geometry.midpoint(s1, s2)
//         // y - y1 = m(x - x1) -> ax + y = c
//         let a = -1 * m
//         let b = 1
//         let c = -1 * (m * midpoint.x - midpoint.y)

//         return {
//             a: a,
//             b: b,
//             c: c,
//         }
//     }

//     static findCircleEvent = (l, m, r, currentY) => {
//         if (!l || !m || !r || (l.p.x === r.p.x && l.p.y === r.p.y)) return null

//         let e1 = Geometry.bisect(l.p, m.p)
//         let e2 = Geometry.bisect(m.p, r.p)

//         // Parallel vertical lines.
//         if (e2.a - e1.a === 0) return null

//         // Intersection of two lines.
//         let x = (e2.c - e1.c) / (e2.a - e1.a)
//         let y = e1.c - e1.a * x

//         let center = { x: x, y: y }

//         // Breakpoint divergence tests.
//         // Has difficulty with degenerate situation where the sites are on a common horizontal line.
//         // See http://davidpritchard.org/graphics/msc_courses/cpsc516/report.pdf
//         if (l.p.y > m.p.y) {
//             if (center.x < l.p.x) return null
//         } else {
//             if (center.x > m.p.x) return null
//         }

//         if (m.p.y > r.p.y) {
//             if (center.x < m.p.x) return null
//         } else {
//             if (center.x > r.p.x) return null
//         }

//         // We are converging!
//         let radius = Geometry.distance(center, m.p) // or l.p or r.p

//         if ((center.y + radius).toFixed(6) <= parseFloat(currentY.toFixed(6))) return null

//         return { x: center.x, y: center.y + radius, r: radius }
//     }

//     /*
//      * lp: the point representing the left region (the point responsible for the arc left of the breakpoint)
//      * rp: the point representing the right region
//      * bisector: perpendicular bisector of lp and rp
//      * z: the y position of the sweep line
//      *
//      * Returns the x position of the breakpoint of interest (either the left or right breakpoint)
//      *
//      * The problem: Given two points A&B and the position of a horizontal line y=z, find the
//      * centers of the two circles that pass through A&B and is tangent w.r.t to line y=z.
//      * If A and B are Voronoi sites, and z is the y position of the sweep line, the centers of
//      * these circles represent the current position of the two breakpoints. Note that the
//      * perpendicular bisector of A&B must pass through the center of both circles.
//      *
//      * Here is an explanation of the mathmetics which has been precalculated for the general case:
//      *  http://math.stackexchange.com/questions/101986/how-to-determine-an-equation-of-a-circle-using-a-line-and-two-points-on-a-circle
//      */
//     static breakpoint = (lp, rp, bisector, z) => {
//         let x = lp.x
//         let y = lp.y
//         let a = bisector.a
//         let c = bisector.c

//         // Solve qA*x^2 + qB*x + qC = 0
//         let qA = -1
//         let qB = 2 * z * a + 2 * x - 2 * y * a
//         let qC = parseFloat((-1 * (x * x + y * y - z * z + 2 * z * c - 2 * y * c)).toFixed(4))
//         let roots = Geometry.roots(qA, qB, qC)

//         if (lp.y > rp.y) return Math.max(roots.r1, roots.r2)
//         else return Math.min(roots.r1, roots.r2)
//     }

//     //Quadratic equation forumla.
//     static roots = (a, b, c) => {
//         let r1 = (b * -1 + Math.sqrt(Math.pow(b, 2) - 4 * a * c)) / (2 * a)
//         let r2 = (b * -1 - Math.sqrt(Math.pow(b, 2) - 4 * a * c)) / (2 * a)
//         return { r1: r1, r2: r2 }
//     }

//     static midpoint = (a, b) => {
//         return {
//             x: Math.abs((a.x + b.x) / 2),
//             y: Math.abs((a.y + b.y) / 2),
//         }
//     }

//     static reciprocal = (a) => {
//         return -1 * (1 / a)
//     }

//     static slope = (a, b) => {
//         return (a.y - b.y) / (a.x - b.x)
//     }

//     static distance = (s, t) => {
//         let dx = s.x - t.x,
//             dy = s.y - t.y
//         return Math.sqrt(dx * dx + dy * dy)
//     }
// }

// export class VoronoiDiagram {
//     constructor(width, height, params) {
//         this.width = width
//         this.height = height
//         this.params = params
//         this.sites = []
//         this.crossed = []

//         this.random = prngd(this.params.seed)
//         for (let i = 0; i < this.params.cells; i++) this.sites.push(this.makeRandomPoint())

//         this.queue = new PriorityQueue((a, b) => a.y > b.y || (a.y == b.y && a.x > b.x))
//         this.beach = new BeachLine()

//         this.sites = this.sites.sort((a, b) => (a.y < b.y ? -1 : a.y > b.y ? 1 : a.x < b.x ? -1 : a.x > b.x ? 1 : 0))
//         for (const p of this.sites) {
//             p.type = 'site'
//             this.queue.insert(p)
//         }

//         while (!this.queue.empty()) {
//             let event = this.queue.pop()
//             if (event.type === 'site') this.handleSiteEvent(event)
//             if (event.type === 'circle') this.handleCircleEvent(event)
//         }

//         this.voronoi = doVoronoi(this.sites.map((v) => [v.x, v.y]))
//         console.log(this.voronoi)
//     }

//     makeRandomPoint = (
//         x = this.random() * this.width,
//         y = this.random() * this.height,
//         v = {
//             x: (this.random() - 0.5) * 2,
//             y: (this.random() - 0.5) * 2,
//         }
//     ) => new Point(x, y, this.width, this.height, v)

//     drawCircle = (ctx, x, y, color) => {
//         ctx.beginPath()
//         ctx.arc(x, y, 4, 0, 2 * Math.PI)
//         ctx.fillStyle = color
//         ctx.fill()

//         ctx.fillStyle = 'black'
//         ctx.fillText('(' + Math.round(x) + ',' + Math.round(y) + ')', x + 3, y)
//     }

//     show(ctx, params) {
//         this.params = params
//         this.sites.forEach((p) => this.drawCircle(ctx, p.x, p.y, 'blue'))
//         this.crossed.forEach((p) => this.drawCircle(ctx, p.x, p.y - p.r, 'red'))
//         this.voronoi.forEach((cell) => {
//             ctx.beginPath()
//             cell.forEach((point, index) => {
//                 if (index === 0) {
//                     ctx.moveTo(point[0], point[1])
//                 } else {
//                     ctx.lineTo(point[0], point[1])
//                 }
//             })
//             ctx.closePath()
//             ctx.stroke()
//         })
//     }

//     attach(canvas) {
//         canvas.addEventListener('mousedown', (e) => {
//             window.restart()
//             const r = canvas.getBoundingClientRect()
//             const p = this.makeRandomPoint(e.clientX - r.left, e.clientY - r.top)
//             this.sites.push(p)
//             this.colors.push(xyToRgb(p.x, p.y, this.width, this.height))
//             window.resume(false)
//         })
//     }

//     //  #######################################################

//     handleSiteEvent(event) {
//         if (this.beach.empty()) {
//             this.beach.initialInsert({ x: event.x, y: event.y })
//             return
//         }

//         let brokenArc = this.beach.search(event)

//         if (brokenArc.event) {
//             queue.remove(brokenArc.event)
//             brokenArc.event = null
//         }

//         let edge = Geometry.bisect(brokenArc.p, event)
//         let arcs = this.beach.split(brokenArc, event, edge)
//         let rr = this.beach.rightArc(arcs.r)
//         let ll = this.beach.leftArc(arcs.l)
//         let r = arcs.r,
//             m = arcs.m,
//             l = arcs.l

//         // We now have 5 consecutive arcs: ll, l, m, r, rr
//         // The triples we must check are (ll, l, m) and (m, r, rr).
//         // We don't need to check (l, m, r) because l and r come from the same parabola.

//         // Check the triple where m is the left arc.
//         this.checkTriple(m, r, rr, event.y)
//         // Check the triple where m is the right arc.
//         this.checkTriple(ll, l, m, event.y)
//     }

//     handleCircleEvent(event) {
//         let deleteArc = event.arc
//         this.crossed.push(event)
//         let neighborsOfRemoved = this.beach.remove(deleteArc, event.y)
//         let l = neighborsOfRemoved.l
//         let r = neighborsOfRemoved.r

//         // Delete circle events on l and r. Note that m was already handled because we pulled it from the queue.
//         if (l.event) this.queue.remove(l.event)
//         l.event = null
//         if (r.event) this.queue.remove(r.event)
//         r.event = null

//         let rr = this.beach.rightArc(r)
//         let ll = this.beach.leftArc(l)

//         // We now have 4 consecutive arcs: ll, l, r, rr (m was removed).
//         // The triples we must check are (ll, l, r) and (l, r, rr).

//         // Check the triple where l is the middle arc.
//         this.checkTriple(ll, l, r, event.y)
//         // Check the triple where r is the middle arc.
//         this.checkTriple(l, r, rr, event.y)
//     }

//     // Checks for a valid circle event among three arcs and adds to the queue.
//     checkTriple(l, m, r, y) {
//         let validCircleEvent = Geometry.findCircleEvent(l, m, r, y)
//         if (validCircleEvent) {
//             validCircleEvent.type = 'circle'
//             validCircleEvent.arc = m
//             let newEvent = this.queue.insert(validCircleEvent)
//             m.event = newEvent
//         }
//     }
// }
