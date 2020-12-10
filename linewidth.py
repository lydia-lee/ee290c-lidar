import numpy as np

def delta_theta(m, lambda_0, lambda_delta, Lambda):
	'''
	Inputs:
		m: Diffraction order.
		lambda_0: Center wavelength (m).
		lambda_delta: Plus/minus wavelength amount to add to and subtract from
			the initial lambda_0 (m).
		Lambda: Phased array reset pitch (m).
	Returns:
		theta_delta: The difference in angle from the nominal value
			that result from lambda_0 +/- lambda_delta.
	'''
	theta_0 = np.arcsin(m*lambda_0/Lambda)
	theta_up = np.arcsin(m*(lambda_0+lambda_delta)/Lambda)
	theta_down = np.arcsin(m*(lambda_0-lambda_delta)/Lambda)

	theta_delta = (theta_0-theta_down, theta_up-theta_0)

	return theta_delta

def delta_f(lambda_0, lambda_delta):
	'''
	Inputs:
		lambda_0: Starting wavelength (m).
		lambda_delta: The change in wavelength, can be positive or negative
			(m).
	Returns:
		f_delta: The corresponding increase or decrease in frequency (Hz)
			corresponding to the change in wavelength.
	'''
	c = 3e8
	f_0 = c/lambda_0
	f_delta = c/(lambda_0 + lambda_delta) - f_0

	return f_delta