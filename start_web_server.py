#!/usr/bin/env python3
"""
AmberFlow - Entry point when running from project root.
Uses the amberflow package. For installed package: use `amberflow` or `python -m amberflow`.
"""

from amberflow.app import app

if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=7860)
