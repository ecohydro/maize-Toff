#! /usr/bin/env python

import unittest
import shutil
import sys

TESTS_COVERAGE_REPORT_DIR = 'cov/'


class CoverageAnalysis(Command):
    """
    Code coverage analysis command.
    """

    description = "run test coverage analysis"
    user_options = [
        ('erase', None, "remove all existing coverage results"),
        ('branch', 'b', 'measure branch coverage in addition to statement coverage'),  # NOQA
        ('test-module=', 't', "explicitly specify a module to test (e.g. 'dendropy.test.test_containers')"),  # NOQA
        ('no-annotate', None, "do not create annotated source code files"),
        ('no-html', None, "do not create HTML report files"),
    ]

    def initialize_options(self):
        """
        Initialize options to default values.
        """
        self.test_module = None
        self.branch = False
        self.erase = False
        self.no_annotate = False
        self.no_html = False
        self.omit_prefixes = ['farm/test']

    def finalize_options(self):
        pass

    def run(self):
        """
        Main command implementation.
        """

        if self.erase:
            try:
                shutil.rmtree(TESTS_COVERAGE_DIR)
            except:
                pass
        else:
            if self.test_module is None:
                test_suite = get_test_suite()
            else:
                test_suite = get_test_suite([self.test_module])
            runner = unittest.TextTestRunner()
            cov = coverage.coverage(branch=self.branch)
            cov.start()
            runner.run(test_suite)
            cov.stop()
            if not self.no_annotate:
                cov.annotate(omit_prefixes=self.omit_prefixes,
                             directory=TESTS_COVERAGE_SOURCE_DIR)
            if not self.no_html:
                cov.html_report(omit_prefixes=self.omit_prefixes,
                                directory=TESTS_COVERAGE_REPORT_DIR)
            cov.report(omit_prefixes=self.omit_prefixes)
