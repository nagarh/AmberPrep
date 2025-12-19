#!/usr/bin/env python3
"""
AmberFlow - Hugging Face Spaces Entry Point
This file serves as the entry point for Hugging Face Spaces deployment.
It imports and runs the existing Flask application from python/app.py
"""

import os
import sys
import subprocess

# Add the python directory to the path so we can import the Flask app
sys.path.append(os.path.join(os.path.dirname(__file__), 'python'))

# Import the Flask app from the existing python/app.py
from python.app import app

if __name__ == "__main__":
    # Run the Flask app
    app.run(debug=False, host='0.0.0.0', port=7860)
