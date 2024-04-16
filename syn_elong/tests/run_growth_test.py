from memote.suite.cli.reports import snapshot


if __name__ == '__main__':

    args = [
        '../new_syn.xml',
        '--filename', 'growth_report_new.html',

        '--experimental', '../data/experiments.yml',
        '--exclusive', 'test_growth',
        '--pytest-args', '--tb=long'


    ]

    snapshot(
       args
    )
