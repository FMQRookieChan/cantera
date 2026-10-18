classdef ctTestNDodecaneCombustor2D < ctTestCase

    methods (TestMethodSetup)

        function addSampleFolder(self)
            thisFile = mfilename('fullpath');
            canteraRoot = fileparts(fileparts(fileparts(thisFile)));
            samplePath = fullfile(canteraRoot, 'samples', 'matlab');
            addpath(samplePath);
        end

    end

    methods (TestMethodTeardown)

        function removeSampleFolder(self)
            thisFile = mfilename('fullpath');
            canteraRoot = fileparts(fileparts(fileparts(thisFile)));
            samplePath = fullfile(canteraRoot, 'samples', 'matlab');
            rmpath(samplePath);
            close all
        end

    end

    methods (Test)

        function testFixedMeshAndFiniteOutputs(self)
            results = ndodecane_combustor_2d(struct(...
                'makePlots', false, ...
                'tEnd', 2.0e-4, ...
                'forceFallbackChemistry', true));

            self.verifyEqual(results.metadata.meshShape, [70, 180]);
            self.verifyEqual(size(results.single.Yfuel), [70, 180]);
            self.verifyEqual(size(results.dual.T), [70, 180]);

            self.verifyTrue(all(all(isfinite(results.single.T))));
            self.verifyTrue(all(all(isfinite(results.single.P))));
            self.verifyTrue(all(all(isfinite(results.dual.T))));
            self.verifyTrue(all(all(isfinite(results.dual.P))));

            self.verifyGreaterThanOrEqual(results.single.mdotFuelBurned, 0.0);
            self.verifyGreaterThanOrEqual(results.dual.mdotFuelBurned, 0.0);
            self.verifyGreaterThan(results.single.outletTemperatureMean, 250.0);
            self.verifyGreaterThan(results.dual.outletPressureMean, 0.0);
        end

        function testSingleDualInjectionSetup(self)
            results = ndodecane_combustor_2d(struct(...
                'makePlots', false, ...
                'tEnd', 2.0e-4, ...
                'forceFallbackChemistry', true));

            self.verifyGreaterThan(max(results.single.Yfuel(:, 1)), 0.0, ...
                'Single-inlet case should inject fuel at x=1 boundary.');

            self.verifyLessThanOrEqual(max(results.dual.Yfuel(:, 1)), 1.0e-12, ...
                'Dual case inlet boundary should not be the primary fuel source.');

            self.verifyGreaterThan(max(max(results.dual.Yfuel(:, 2:end))), 0.0, ...
                'Dual-inlet case should inject fuel into interior configured regions.');
        end

        function testMechanismPathResolution(self)
            results = ndodecane_combustor_2d(struct(...
                'makePlots', false, ...
                'tEnd', 1.0e-4, ...
                'forceFallbackChemistry', true));

            mechanismPath = results.metadata.mechanismPath;
            [~, mechName, mechExt] = fileparts(mechanismPath);

            self.verifyTrue(isfile(mechanismPath));
            self.verifyEqual([mechName, mechExt], 'nDodecane_Reitz.yaml');
            self.verifyEqual(results.metadata.indexing, '[radial, axial]');
        end

    end

end
