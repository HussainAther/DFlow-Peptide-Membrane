class AmphiphileLimiter:
    def __init__(self, initial_pool=500):
        self.pool = initial_pool
        self.used = 0

    def can_grow(self, peptides_needed):
        return (self.pool - self.used) >= peptides_needed

    def consume(self, amount):
        self.used += amount
        self.pool = max(0, self.pool - amount)

    def reset(self):
        self.used = 0

