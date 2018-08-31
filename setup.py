from setuptools import setup

setup(
    name="dparser",
    version='0.1',
    py_modules=['DataParser'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        dparser=DataParser:cli
    ''',
)