# Contributing to Sprint Satellite

Thank you for your interest in contributing to the Sprint Satellite project! This document provides guidelines and information for contributors.

## Project Goals

The Sprint Satellite project aims to:
1. Demonstrate the viability of high-power, short-duration CubeSat missions
2. Provide accurate physics simulations of electrodynamic tether systems
3. Compare MHD satellite performance against traditional solar-powered satellites
4. Enable research and development of novel satellite power generation concepts

## How to Contribute

### Reporting Issues

Before creating an issue, please:
1. Check if the issue has already been reported
2. Search the existing issues for similar problems
3. Provide detailed information including:
   - Python version
   - Operating system
   - Error messages
   - Steps to reproduce

### Suggesting Enhancements

We welcome suggestions for:
- Additional satellite configurations
- Improved physics models
- Enhanced visualization options
- Performance optimizations
- Documentation improvements

### Code Contributions

#### Development Setup

1. Fork the repository
2. Clone your fork locally
3. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
4. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

#### Coding Standards

- Follow PEP 8 style guidelines
- Use meaningful variable and function names
- Add docstrings to all functions and classes
- Include type hints where appropriate
- Write unit tests for new functionality

#### Testing

Before submitting a pull request:
1. Run the simulation to ensure it works:
   ```bash
   python main.py
   ```
2. Verify that plots are generated correctly
3. Check that the summary CSV is created
4. Test with different parameter configurations

#### Pull Request Process

1. Create a feature branch from `main`
2. Make your changes with clear commit messages
3. Update documentation if needed
4. Test your changes thoroughly
5. Submit a pull request with:
   - Clear description of changes
   - Reference to any related issues
   - Screenshots of new visualizations (if applicable)

## Areas for Contribution

### Physics Models
- More sophisticated atmospheric models
- Advanced magnetic field calculations
- Improved eclipse detection algorithms
- Additional perturbation forces

### Visualization
- Interactive plots
- 3D orbital visualization
- Real-time simulation display
- Export to different formats

### Performance
- Parallel processing for multiple satellites
- Optimized API calls
- Caching improvements
- Memory usage optimization

### Documentation
- API documentation
- Tutorial notebooks
- Video demonstrations
- Academic paper references

## Code of Conduct

- Be respectful and inclusive
- Focus on constructive feedback
- Help others learn and contribute
- Follow the project's technical decisions

## Getting Help

- Open an issue for bugs or questions
- Join discussions in existing issues
- Check the README for basic usage
- Review the code comments for implementation details

## License

By contributing to this project, you agree that your contributions will be licensed under the MIT License.

Thank you for contributing to the Sprint Satellite project! 