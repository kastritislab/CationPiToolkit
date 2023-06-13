from setuptools import setup, find_packages

setup(
    name='CationPiToolkit',
    version='0.1.0',
    packages=find_packages(),
    url='https://github.com/ctueting/CationPiToolkit', 
    license='MIT',
    author='Christian Tüting',
    author_email='christian.tueting@biochemtech.uni-halle.de',
    description='A toolkit for analysing a protein structure for potential cation pi interactions.',
    long_description=open('README.md').read(), 
    long_description_content_type="text/markdown", 
    install_requires=[
        'numpy',
        'pandas',
        'tqdm', 
    ],
    python_requires='>=3.6', 
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'CationPiToolkit=CationPiToolkit.CationPiToolkit:main',
        ],
    },
)
