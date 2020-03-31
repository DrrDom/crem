import setuptools
from os import path
import crem

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name="crem",
    version=crem.__version__,
    author="Pavel Polishchuk",
    author_email="pavel_polishchuk@ukr.net",
    description="CReM: chemically reasonable mutations framework",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DrrDom/crem",
    packages=['crem'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    python_requires='>=3.6',
    extras_require={
        'rdkit': ['rdkit>=2017.09'],
    },
    entry_points={'console_scripts':
                      ['fragmentation = crem.fragmentation:entry_point',
                       'frag_to_env = crem.frag_to_env_mp:entry_point',
                       'env_to_db = crem.import_env_to_db:entry_point',
                       'guacamol_test = crem.guacamol_crem_test:entry_point']},
    scripts=['crem/scripts/crem_create_frag_db.sh']
)
