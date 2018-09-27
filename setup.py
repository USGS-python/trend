import setuptools

with open("README.md","r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='trend',
    version='0.1',
    description='Algorithms for detecting trends in time-series data.',
    url='https://github.com/USGS-python/trend',
    author='Timothy Hodson',
    author_email='thodson@usgs.gov',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    packages=setuptools.find_packages(),
    zip_safe=False
)
