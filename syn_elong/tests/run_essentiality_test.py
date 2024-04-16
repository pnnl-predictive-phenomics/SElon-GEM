from memote.suite.cli.reports import snapshot


if __name__ == '__main__':

    args = [
        '../syn_elong.xml',
        '--filename', 'essential_only_report.html',

        # '../new_syn.xml',
        # '--filename', 'essential_only_report2.html',


        '--experimental', '../data/experiments.yml',
        # '--exclusive', 'test_essentiality',
        '--pytest-args', '--tb=long'


    ]

    snapshot(
       args
    )
