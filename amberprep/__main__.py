"""Run the AmberPrep web server. Use: python -m amberprep or amberprep"""

from amberprep.app import app


def main():
    app.run(debug=False, host="0.0.0.0", port=7860)


if __name__ == "__main__":
    main()
