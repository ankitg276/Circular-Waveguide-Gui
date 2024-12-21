# Circular-Waveguide-Gui
classdef TE_Waveguide_App < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        CutoffFrequencyLabel            matlab.ui.control.Label
        ModeNumbermEditField            matlab.ui.control.NumericEditField
        ModeNumbermLabel                matlab.ui.control.Label
        RadiusofthewaveguiderEditField  matlab.ui.control.NumericEditField
        RadiusofthewaveguiderLabel      matlab.ui.control.Label
        CalculateButton                 matlab.ui.control.Button
        ModeNumbernEditField            matlab.ui.control.NumericEditField
        ModeNumbernLabel                matlab.ui.control.Label
        RelativePermeabilitymu_rEditField  matlab.ui.control.NumericEditField
        RelativePermeabilitymu_rLabel   matlab.ui.control.Label
        RelativePermittivityepsilon_rEditField  matlab.ui.control.NumericEditField
        RelativePermittivityepsilon_rLabel  matlab.ui.control.Label
        UIAxes4                         matlab.ui.control.UIAxes
        UIAxes3                         matlab.ui.control.UIAxes
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes                          matlab.ui.control.UIAxes
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(app, event)
  % Get input values from the UI
            r = app.RadiusofthewaveguiderEditField.Value;
            epsilon_r = app.RelativePermittivityepsilon_rEditField.Value;
            mu_r = app.RelativePermeabilitymu_rEditField.Value;
            m = app.ModeNumbermEditField.Value;
            n = app.ModeNumbernEditField.Value;

            % Input validation
            if r <= 0
                uialert(app.UIFigure, 'Radius must be a positive value.', 'Input Error');
                return;
            end
            if epsilon_r <= 0
                uialert(app.UIFigure, 'Relative permittivity must be a positive value.', 'Input Error');
                return;
            end
            if mu_r <= 0
                uialert(app.UIFigure, 'Relative permeability must be a positive value.', 'Input Error');
                return;
            end
            if m <= 0 || m ~= floor(m)
                uialert(app.UIFigure, 'Mode number m must be a positive integer.', 'Input Error');
                return;
            end
            if n <= 0 || n ~= floor(n)
                uialert(app.UIFigure, 'Mode number n must be a positive integer.', 'Input Error');
                return;
            end

            % Constants
            c = 3e8; % Speed of light in vacuum (m/s)
            epsilon_0 = 8.854e-12; % Permittivity of free space (F/m)
            mu_0 = 4 * pi * 1e-7; % Permeability of free space (H/m)

            % Derived constants
            epsilon = epsilon_r * epsilon_0;
            mu = mu_r * mu_0;

            % Calculate the nth root of the m-th order Bessel function of the first kind
            bessel_func = @(x) besselj(m, x); % Bessel function J_m(x)
            chi_mn = zeros(1, n); % Preallocate for n roots

            % Find the first n roots of the Bessel function for the given m
            root_guess = m; % Initial guess (roots roughly spaced by pi)
            for k = 1:n
                try
                    chi_mn(k) = fzero(bessel_func, [root_guess, root_guess + pi]);
                catch
                    uialert(app.UIFigure, 'Unable to find Bessel function roots. Try adjusting your parameters.', 'Calculation Error');
                    app.LoadingLabel.Visible = 'off';
                    return;
                end
                root_guess = chi_mn(k) + pi; % Update guess for the next root
            end
            chi_mn = chi_mn(n); % Use the nth root

            % Calculate the cutoff frequency
            f_c = (c / (2 * pi * r * sqrt(epsilon_r * mu_r))) * chi_mn;

            % Display the cutoff frequency on the UI
            app.CutoffFrequencyLabel.Text = ['Cutoff Frequency: ' num2str(f_c / 1e9, '%.2f') ' GHz'];

            % Dynamically adjust grid resolution based on waveguide radius
            max_resolution = 500;
            min_resolution = 50; % Minimum resolution for performance
            resolution = round(300 * (r / 10)); % Adjust resolution based on radius
            resolution = max(min_resolution, min(resolution, max_resolution)); % Limit resolution within a range

            % Define the 2D grid for x-y plane
            x = linspace(-r, r, resolution);
            y = linspace(-r, r, resolution);
            [X, Y] = meshgrid(x, y);
            [theta, rho] = cart2pol(X, Y);

            % Compute the cutoff wavenumber
            kc = chi_mn / r;

            % Angular frequency for this cutoff frequency
            omega = 2 * pi * f_c;

            % Calculate field components for TE_mn mode
            Hz = besselj(m, kc * rho) .* cos(m * theta);

            % Electric field components in cylindrical coordinates
            [Hz_rho, ~] = gradient(Hz, x(2)-x(1));
            Er = -kc ./ (omega * epsilon) .* Hz_rho;
            Ephi = (m ./ (rho * omega * epsilon)) .* Hz;

            % Convert cylindrical to Cartesian components for plotting in x-y plane
            Ex = Er .* cos(theta) - Ephi .* sin(theta);
            Ey = Er .* sin(theta) + Ephi .* cos(theta);

            % Normalize field vectors for clearer visualization
            magnitude = max(sqrt(Ex.^2 + Ey.^2), [], 'all');
            Ex = Ex / magnitude;
            Ey = Ey / magnitude;

            % Plot electric field lines in the x-y plane
            axes(app.UIAxes);
            quiver(app.UIAxes, X, Y, Ex, Ey, 2); % Increased scaling factor for clarity
            axis(app.UIAxes, 'equal');
            title(app.UIAxes, ['Electric Field Lines in x-y Plane for TE_{' num2str(m) num2str(n) '} Mode']);
            xlabel(app.UIAxes, 'x (m)');
            ylabel(app.UIAxes, 'y (m)');
            legend(app.UIAxes, 'Electric Field Lines', 'Location', 'Best');

            % Enable zoom and pan interaction for the axes
            zoom(app.UIAxes, 'on');
            pan(app.UIAxes, 'on');

            % Magnetic field components in cylindrical coordinates
            [Er_x, ~] = gradient(Er, x(2)-x(1));
            Hr = -kc ./ (omega * mu) .* Er_x;
            Hphi = (m ./ (rho * omega * mu)) .* Er;

            % Convert cylindrical to Cartesian components for plotting in x-y plane
            Hx = Hr .* cos(theta) - Hphi .* sin(theta);
            Hy = Hr .* sin(theta) + Hphi .* cos(theta);

            % Normalize magnetic field vectors
            magnitude_h = max(sqrt(Hx.^2 + Hy.^2), [], 'all');
            Hx = Hx / magnitude_h;
            Hy = Hy / magnitude_h;

            % Plot magnetic field lines in the x-y plane
            axes(app.UIAxes3);
            quiver(app.UIAxes3, X, Y, Hx, Hy, 2);
            axis(app.UIAxes3, 'equal');
            title(app.UIAxes3, ['Magnetic Field Lines in x-y Plane for TE_{' num2str(m) num2str(n) '} Mode']);
            xlabel(app.UIAxes3, 'x (m)');
            ylabel(app.UIAxes3, 'y (m)');
            legend(app.UIAxes3, 'Magnetic Field Lines', 'Location', 'Best');

            % Enable zoom and pan interaction for the axes
            zoom(app.UIAxes3, 'on');
            pan(app.UIAxes3, 'on');

            % Define the grid for y-z plane
            z = linspace(0, 0.2, resolution);  % Increased range for better visibility
            [Y, Z] = meshgrid(y, z);

            % Calculate and plot electric field in y-z plane
            axes(app.UIAxes2);
            quiver(app.UIAxes2, Y, Z, Ex, Ey, 2);
            axis(app.UIAxes2, 'equal');
            title(app.UIAxes2, 'Electric Field in y-z Plane');
            xlabel(app.UIAxes2, 'y (m)');
            ylabel(app.UIAxes2, 'z (m)');

            % Enable zoom and pan interaction for the axes
            zoom(app.UIAxes2, 'on');
            pan(app.UIAxes2, 'on');

            % Magnetic field components for the y-z plane
            [Er_y, ~] = gradient(Er, y(2)-y(1));  % Gradient along y
            Hz_y = -kc ./ (omega * mu) .* Er_y;
            Hphi_y = (m ./ (y * omega * mu)) .* Er;  % Modify this accordingly for y-z plane

            % Convert cylindrical to Cartesian coordinates for magnetic field
            Hx_y = Hz_y .* cos(theta) - Hphi_y .* sin(theta);
            Hy_y = Hz_y .* sin(theta) + Hphi_y .* cos(theta);

            % Normalize magnetic field components
            magnitude_h_y = max(sqrt(Hx_y.^2 + Hy_y.^2), [], 'all');
            Hx_y = Hx_y / magnitude_h_y;
            Hy_y = Hy_y / magnitude_h_y;

            % Plot magnetic field in y-z plane
            axes(app.UIAxes4);
            quiver(app.UIAxes4, Y, Z, Hx_y, Hy_y, 2);
            axis(app.UIAxes4, 'equal');
            title(app.UIAxes4, 'Magnetic Field Lines in y-z Plane');
            xlabel(app.UIAxes4, 'y (m)');
            ylabel(app.UIAxes4, 'z (m)');

            % Enable zoom and pan interaction for the axes
            zoom(app.UIAxes4, 'on');
            pan(app.UIAxes4, 'on');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Electric Field Lines for TEmn')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [317 207 300 185];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, 'Electric Field Lines for TEmn')
            xlabel(app.UIAxes2, 'Y')
            ylabel(app.UIAxes2, 'Z')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.Position = [317 11 300 185];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.UIFigure);
            title(app.UIAxes3, 'Magnetic Field Lines for TEmn')
            xlabel(app.UIAxes3, 'X')
            ylabel(app.UIAxes3, 'Y')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.Position = [18 207 300 185];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.UIFigure);
            title(app.UIAxes4, 'Magnetic Field Lines for TEmn')
            xlabel(app.UIAxes4, 'Y')
            ylabel(app.UIAxes4, 'Z')
            zlabel(app.UIAxes4, 'Z')
            app.UIAxes4.Position = [9 11 300 185];

            % Create RelativePermittivityepsilon_rLabel
            app.RelativePermittivityepsilon_rLabel = uilabel(app.UIFigure);
            app.RelativePermittivityepsilon_rLabel.HorizontalAlignment = 'right';
            app.RelativePermittivityepsilon_rLabel.Position = [201 449 130 22];
            app.RelativePermittivityepsilon_rLabel.Text = 'Relative Permittivity (epsilon_r)';

            % Create RelativePermittivityepsilon_rEditField
            app.RelativePermittivityepsilon_rEditField = uieditfield(app.UIFigure, 'numeric');
            app.RelativePermittivityepsilon_rEditField.Position = [330 449 78 22];

            % Create RelativePermeabilitymu_rLabel
            app.RelativePermeabilitymu_rLabel = uilabel(app.UIFigure);
            app.RelativePermeabilitymu_rLabel.HorizontalAlignment = 'right';
            app.RelativePermeabilitymu_rLabel.Position = [422 449 126 22];
            app.RelativePermeabilitymu_rLabel.Text = 'Relative Permeability (mu_r)';

            % Create RelativePermeabilitymu_rEditField
            app.RelativePermeabilitymu_rEditField = uieditfield(app.UIFigure, 'numeric');
            app.RelativePermeabilitymu_rEditField.Position = [557 449 72 22];

            % Create ModeNumbernLabel
            app.ModeNumbernLabel = uilabel(app.UIFigure);
            app.ModeNumbernLabel.HorizontalAlignment = 'right';
            app.ModeNumbernLabel.Position = [169 417 104 22];
            app.ModeNumbernLabel.Text = 'Mode Number (n)';

            % Create ModeNumbernEditField
            app.ModeNumbernEditField = uieditfield(app.UIFigure, 'numeric');
            app.ModeNumbernEditField.Position = [280 417 62 22];

            % Create CalculateButton
            app.CalculateButton = uibutton(app.UIFigure, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.Position = [350 417 73 23];
            app.CalculateButton.Text = 'Calculate';

            % Create RadiusofthewaveguiderLabel
            app.RadiusofthewaveguiderLabel = uilabel(app.UIFigure);
            app.RadiusofthewaveguiderLabel.HorizontalAlignment = 'right';
            app.RadiusofthewaveguiderLabel.Position = [6 449 106 22];
            app.RadiusofthewaveguiderLabel.Text = 'Radius of the waveguide (r)';

            % Create RadiusofthewaveguiderEditField
            app.RadiusofthewaveguiderEditField = uieditfield(app.UIFigure, 'numeric');
            app.RadiusofthewaveguiderEditField.Position = [111 449 72 22];

            % Create ModeNumbermLabel
            app.ModeNumbermLabel = uilabel(app.UIFigure);
            app.ModeNumbermLabel.HorizontalAlignment = 'right';
            app.ModeNumbermLabel.Position = [1 417 106 22];
            app.ModeNumbermLabel.Text = 'Mode Number (m)';

            % Create ModeNumbermEditField
            app.ModeNumbermEditField = uieditfield(app.UIFigure, 'numeric');
            app.ModeNumbermEditField.Position = [111 417 59 22];

            % Create CutoffFrequencyLabel
            app.CutoffFrequencyLabel = uilabel(app.UIFigure);
            app.CutoffFrequencyLabel.Position = [436 417 193 22];
            app.CutoffFrequencyLabel.Text = 'Cutoff Frequency';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = TE_Waveguide_App

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
