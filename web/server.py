from flask import *
app = Flask(__name__)


# Note: need to annotante number of loci (MSAT, SNP)


def config():
    pass


@app.route('/models')
def show_models():
    pass


@app.route('/model')
def do_model():
    pass


@app.route('/simulations')
def show_simulations():
    pass


@app.route('/analysis')
def show_analysis():
    pass


@app.route('/')
def index():
        return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
