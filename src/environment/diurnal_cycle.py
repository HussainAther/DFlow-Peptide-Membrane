class DiurnalCycle:
    def __init__(self, uv_max=1.0, hydration_peak=1.0, period=2):
        self.phase = 0
        self.period = period  # in "steps" (e.g., days)
        self.uv_max = uv_max
        self.hydration_peak = hydration_peak

    def step(self):
        self.phase = (self.phase + 1) % self.period

    def get_uv_level(self):
        return self.uv_max if self.phase == 0 else 0.1

    def get_hydration_level(self):
        return self.hydration_peak if self.phase == 1 else 0.2

    def get_temperature(self):
        return 25 if self.phase == 0 else 15  # simple model

