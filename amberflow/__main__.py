"""Run the AmberFlow web server. Use: python -m amberflow or amberflow"""

from amberflow.app import app


def main():
    app.run(debug=False, host="0.0.0.0", port=7860)


if __name__ == "__main__":
    main()
