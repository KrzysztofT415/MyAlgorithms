class PriorityQueueData {
    constructor(data, priority) {
        this.data_stored = data
        this.priority = priority
    }
}

export class PriorityQueue {
    constructor() {
        this.data_stored = []
    }

    get isEmpty() {
        return this.data_stored.length === 0
    }

    get front() {
        return this.data_stored[0].data_stored
    }

    insert = (data, priority) => {
        const id = this.data_stored.findIndex((cmp) => cmp.priority > priority)
        const node = new PriorityQueueData(data, priority)
        if (id !== -1) this.data_stored.splice(id, 0, node)
        else this.data_stored.push(node)
    }

    remove = (data) => {
        const id = this.data_stored.findIndex((cmp) => cmp.data_stored == data)
        if (id !== -1) this.data_stored.splice(id, 1)
    }

    update = (data, new_priority) => {
        this.remove(data)
        this.insert(data, new_priority)
    }

    pop = () => this.data_stored.shift().data_stored
}
