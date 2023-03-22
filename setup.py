from distutils.core import setup


def main():
    setup(name='syn_elong',
          version='1.0',
          description='synechococcus_elongatus_pcc_7942',
          author='?',
          author_email='?',
          packages=['syn_elong'],
          requires=['cobra'],
          keywords=['systems', 'biology', 'model', 'rules'],
          classifiers=[
            'Intended Audience :: Science/Research',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            ],
          )


if __name__ == '__main__':
    main()