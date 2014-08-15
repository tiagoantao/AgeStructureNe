import sys
sys.path.append('..')

import os
import sqlite3

from flask import Flask, g, render_template
app = Flask(__name__)

import myUtils

# Load default config and override config from an environment variable
app.config.update(dict(
    DATABASE=os.path.join(app.root_path, 'age.db'),
    DEBUG=True,
))


def connect_db():
    rv = sqlite3.connect(app.config['DATABASE'])
    rv.row_factory = sqlite3.Row
    return rv


def init_db():
    with app.app_context():
        db = get_db()
        with app.open_resource('schema.sql', mode='r') as f:
            db.cursor().executescript(f.read())
        db.commit()


def get_db():
    if not hasattr(g, 'db'):
        g.sqlite_db = connect_db()
    return g.sqlite_db


@app.teardown_appcontext
def close_db(error):
    if hasattr(g, '_db'):
        g.sqlite_db.close()


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
