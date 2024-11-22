"""
Working Test script companion to observation_data.py
"""
    """A Simple Target class"""

    def __init__(self, name, ra, dec):
        self.name = name
        self.ra = ra
        self.dec = dec
        # Can add other properties of Target from database or apis

    def __repr__(self):
        return "Target('{}', '{}', {})".format(self.name, self.ra, self.dec)
"""
