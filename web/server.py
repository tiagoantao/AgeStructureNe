import sys
sys.path.append('..')

from flask import *
app = Flask(__name__)

import myUtils


# Note: need to annotante number of loci (MSAT, SNP)


def config():
    pass


@app.route('/models')
def show_models():
    return render_template('models.html', meta=myUtils.meta)


@app.route('/model/<key>')
def do_model(key):
    return key


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
