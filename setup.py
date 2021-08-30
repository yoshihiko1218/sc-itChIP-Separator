from setuptools import setup, find_packages
  
with open('requirements.txt') as f:
    requirements = f.readlines()
  
long_description = 'Sample Package made for a demo \
      of its making for the GeeksforGeeks Article.'
  
setup(
        name ='scitChIP_Sep',
        version ='0.0.1',
        author="Jiayan(Yoshii) Ma",
        author_email="jim095@ucsd.edu",
        url ='https://github.com/yoshihiko1218/sc-itChIP-Separator',
        description="A package to make file preparation for sc-itChIP sequence files",
        long_description = long_description,
        long_description_content_type ="text/markdown",
        license ='MIT',
        packages = find_packages(),
        entry_points ={
            'console_scripts': [
                'scitChipSeparate = scitChIP_Sep.scitChipSeparate:main'
            ]
        },
        classifiers =[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        keywords ='scitChIP preprocess separator package',
        install_requires = requirements,
        zip_safe = False
)