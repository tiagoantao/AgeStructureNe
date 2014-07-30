import flask
app = flask.Flask(__name__)


# Note: need to annotante number of loci (MSAT, SNP)


@app.route('/')
def index():
    return flask.url_for('models')


@app.route('/config')
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

if __name__ == '__main__':
    app.run()
