function tests = test_v22_fields
% TEST_V22_FIELDS The optional instance fields of PulSeg 2.2.0 in validate_ir and write_seq.
%
% Requires the MATLAB Pulseq toolbox on the path.

tests = functiontests(localfunctions);
end


function setupOnce(testCase)
    testCase.assumeTrue( ...
        exist('mr.Sequence', 'class') == 8 || exist('mr.Sequence', 'file') == 2, ...
        'Pulseq MATLAB toolbox was not found on the MATLAB path.');
    tmpDir = tempname;
    mkdir(tmpDir);
    testCase.TestData.tmpDir = tmpDir;
    seqFile = fullfile(tmpDir, 'demo.seq');
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'demo'));
    create_demo_sequence(seqFile);
    testCase.TestData.ir = pulseg.import(seqFile, 'verbose', false);
end


function teardownOnce(testCase)
    if isfolder(testCase.TestData.tmpDir)
        rmdir(testCase.TestData.tmpDir, 's');
    end
end


function testImportWritesVersion220(testCase)
    testCase.verifyEqual(testCase.TestData.ir.pulseg_version, '2.2.0');
    pulseg.validate_ir(testCase.TestData.ir);
end


function testFlagsNeedVersion220(testCase)
    ir = testCase.TestData.ir;
    ir.pulseg_version = '2.1-alpha';
    ir.execution_stream(1).norot = 1;
    testCase.verifyError(@() pulseg.validate_ir(ir), ?MException);
end


function testFlagMustBeZeroOrOne(testCase)
    ir = testCase.TestData.ir;
    ir.execution_stream(1).nopos = 2;
    testCase.verifyError(@() pulseg.validate_ir(ir), ?MException);
end


function testNorotIsWrittenAsALabelWhereTheInstanceSetsIt(testCase)
    ir = testCase.TestData.ir;
    for k = 1:numel(ir.execution_stream)
        ir.execution_stream(k).norot = double(k == numel(ir.execution_stream));
    end
    seq = pulseg.write_seq(ir, fullfile(testCase.TestData.tmpDir, 'norot.seq'));
    labels = seq.evalLabels('evolution', 'blocks');
    testCase.verifyTrue(isfield(labels, 'NOROT'));
    testCase.verifyEqual(labels.NOROT(1), 0);
    testCase.verifyEqual(labels.NOROT(end), 1);
end


function testRfShimIsWrittenWithEachRfBlock(testCase)
    ir = testCase.TestData.ir;
    for k = 1:numel(ir.execution_stream)
        n_rf = numel(ir.execution_stream(k).rf_amplitude);
        ir.execution_stream(k).rf_shim = repmat([1, 0.5i], n_rf, 1);
    end
    pulseg.validate_ir(ir);
    seq = pulseg.write_seq(ir, fullfile(testCase.TestData.tmpDir, 'shim.seq'));
    testCase.verifyGreaterThan(seq.duration(), 0);
end


function testReadoutPlayedButNotStoredIsRefused(testCase)
    ir = testCase.TestData.ir;
    for k = 1:numel(ir.execution_stream)
        n_adc = numel(ir.execution_stream(k).adc_phase_offset);
        if n_adc > 0
            ir.execution_stream(k).adc_acquire = zeros(n_adc, 1);
            break;
        end
    end
    testCase.verifyError(@() pulseg.write_seq(ir, fullfile(testCase.TestData.tmpDir, 'adc.seq')), ?MException);
end
