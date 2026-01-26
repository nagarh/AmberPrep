#!/usr/bin/env python3
"""
AmberPrep - Entry point when running from project root.
Uses the amberprep package. For installed package: use `amberprep` or `python -m amberprep`.
"""

from amberprep.app import app

if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=7860)
