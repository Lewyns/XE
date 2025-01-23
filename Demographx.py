import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import messagebox
import math

class GraphingCalculator:
    def __init__(self, root):
        self.root = root
        self.root.title('DEMO_GRAPHX')

        # Initialize symbols
        self.x = sp.Symbol('x')
        self.t = sp.Symbol('t')
        self.parameters = {'t': 0}

        # Create widgets
        self.create_widgets()
        
        self.figure = None
        self.ax = None

    def create_widgets(self):
        # Input for mathematical expressions
        tk.Label(self.root, text="Enter mathematical expressions (comma separated):").pack(pady=5)
        self.entry_formula = tk.Entry(self.root, width=50)
        self.entry_formula.pack(pady=5)

        # Input for parameters
        tk.Label(self.root, text="Enter parameters (name=value, comma separated):").pack(pady=5)
        self.entry_parameters = tk.Entry(self.root, width=50)
        self.entry_parameters.pack(pady=5)

        # Calculate button
        calc_button = tk.Button(self.root, text="Calculate", command=self.calculate_result)
        calc_button.pack(pady=5)

        # Display result
        self.label_result = tk.Label(self.root, text="Result: ")
        self.label_result.pack(pady=5)

        # Plot button
        plot_button = tk.Button(self.root, text="Plot Graph", command=self.plot_graph)
        plot_button.pack(pady=5)

        # Zoom In button
        zoom_in_button = tk.Button(self.root, text="Zoom In", command=self.zoom_in)
        zoom_in_button.pack(pady=5, side=tk.LEFT, padx=5)

        # Zoom Out button
        zoom_out_button = tk.Button(self.root, text="Zoom Out", command=self.zoom_out)
        zoom_out_button.pack(pady=5, side=tk.LEFT, padx=5)

        # Reset View button
        reset_button = tk.Button(self.root, text="Reset View", command=self.reset_view)
        reset_button.pack(pady=5, side=tk.LEFT, padx=5)

        # Coordinates label
        self.label_coord = tk.Label(self.root, text="Current coordinates: (x, y)")
        self.label_coord.pack(pady=5)

        # Plot frame
        self.frame_plot = tk.Frame(self.root)
        self.frame_plot.pack(fill=tk.BOTH, expand=True)

    def calculate_result(self):
        try:
            # Parse the parameters
            param_strs = self.entry_parameters.get().split(',')
            self.parameters = {'t': 0}
            for param_str in param_strs:
                if '=' in param_str:
                    name, value = param_str.split('=')
                    self.parameters[name.strip()] = float(value.strip())

            expressions = self.entry_formula.get().split(',')
            results = []
            for expr_str in expressions:
                expr_str = expr_str.replace('y =', '').strip()  # Remove 'y =' if exists

                # Handle expressions like 'y*5 = x**3' or inequalities
                if '=' in expr_str:
                    lhs, rhs = expr_str.split('=')
                    expr = sp.sympify(f"{lhs.strip()} - ({rhs.strip()})")
                    solutions = sp.solve(expr, self.x)
                    results.append(f"Solutions for {expr_str}: {solutions}")
                elif '>' in expr_str or '<' in expr_str:
                    expr = sp.sympify(expr_str)
                    solutions = sp.solve_univariate_inequality(expr, self.x, relational=False)
                    results.append(f"Solutions for {expr_str}: {solutions}")
                else:
                    expr = sp.sympify(expr_str)
                    result = expr.evalf(subs={**self.parameters})
                    results.append(f"{expr_str}: {result}")
            
            self.label_result.config(text="Result: " + ", ".join(results))
        except Exception as e:
            messagebox.showerror("Error", f"Invalid expression or value: {e}")

    def plot_graph(self):
        try:
            # Parse the parameters
            param_strs = self.entry_parameters.get().split(',')
            self.parameters = {'t': 0}
            for param_str in param_strs:
                if '=' in param_str:
                    name, value = param_str.split('=')
                    self.parameters[name.strip()] = float(value.strip())

            # Increase the range of x_vals to [-10π, 10π] for better visualization of trigonometric functions
            x_vals = np.linspace(-10 * np.pi, 10 * np.pi, 1000)

            if self.figure is not None:
                self.figure.clear()

            self.figure, self.ax = plt.subplots()

            expressions = self.entry_formula.get().split(',')
            for expr_str in expressions:
                expr_str = expr_str.replace('y =', '').strip()  # Remove 'y =' if exists

                if '=' in expr_str:  # Handle equations like 'y**2 = x + 5'
                    lhs, rhs = expr_str.split('=')
                    expr = sp.sympify(f"{lhs.strip()} - ({rhs.strip()})")
                    y_solutions = sp.solve(expr, sp.Symbol('y'))
                    for y_expr in y_solutions:
                        f = sp.lambdify(self.x, y_expr, 'numpy')
                        y_vals = f(x_vals)
                        self.ax.plot(x_vals, y_vals, label=str(y_expr))

                elif '>' in expr_str or '<' in expr_str:  # Handle inequalities
                    expr = sp.sympify(expr_str)
                    solutions = sp.solve_univariate_inequality(expr, self.x, relational=False)

                    if isinstance(solutions, sp.Interval):
                        # Convert interval bounds to floats
                        x_min = float(solutions.start)
                        x_max = float(solutions.end)
                        mask = (x_vals >= x_min) & (x_vals <= x_max)
                        self.ax.fill_between(x_vals, -10, 10, where=mask, alpha=0.3)

                    else:
                        # If not an interval, handle other types appropriately
                        numeric_solutions = [float(sol.evalf()) for sol in solutions]
                        for val in numeric_solutions:
                            self.ax.axvline(val, color='r', linestyle='--', label=f'x = {val}')

                else:  # Handle simple expressions like 'y = x**2'
                    expr = sp.sympify(expr_str)
                    expr = expr.subs(sp.Abs, lambda x: math.fabs(x))
                    f = sp.lambdify(self.x, expr, 'numpy')
                    y_vals = f(x_vals)
                    if 'Abs' in str(expr):
                # หาจุดที่ทำให้ค่าภายในสัมบูรณ์เป็น 0
                        inner_expr = expr.args[0]  # สมมติว่าค่าสัมบูรณ์อยู่ที่อาร์กิวเมนต์แรก
                        roots = sp.solve(inner_expr, self.x)
                        roots = [float(root) for root in roots if isinstance(root, sp.Number)]
                # เพิ่มจุดเหล่านี้ลงใน x_vals เพื่อให้กราฟมีความต่อเนื่อง
                        x_vals = np.sort(np.concatenate([x_vals, roots]))
                        y_vals = f(x_vals)
                    intersection_x = sp.solve(expr, self.x)
                    # Handle potential for multiple solutions
                    for x_val in intersection_x:
                        if isinstance(x_val, sp.Symbol):
                            continue  # Skip symbolic solutions
                    y_val = f(x_val)
                    self.ax.scatter(x_val, y_val, marker='o', c='red', label=f'Intersection (x, y) = ({x_val:.2f}, {y_val:.2f})')

                    # Handle discontinuities in functions like tan(x)
                    if 'tan' in expr_str:
                        y_vals = np.ma.masked_where(np.abs(y_vals) > 10, y_vals)

                    self.ax.plot(x_vals, y_vals, label=str(expr))

            self.ax.set_title('Graph of Expressions')
            self.ax.set_xlabel('x')
            self.ax.set_ylabel('y')
            self.ax.legend()
            self.ax.grid(True)

            # Add interactive click events
            self.ax.figure.canvas.mpl_connect('button_press_event', self.on_click)

            # Render the plot
            for widget in self.frame_plot.winfo_children():
                widget.destroy()

            canvas = FigureCanvasTkAgg(self.figure, master=self.frame_plot)
            self.ax.set_xlim([-10, 10])
            self.ax.set_ylim([-10, 10])
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        except Exception as e:
            messagebox.showerror("Error", f"Invalid expression: {e}")

    def on_click(self, event):
        if self.ax is not None and event.inaxes == self.ax:
            xdata, ydata = event.xdata, event.ydata
            self.label_coord.config(text=f"Current coordinates: ({xdata:.2f}, {ydata:.2f})")

    def zoom_in(self):
        if self.ax:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self.ax.set_xlim([x * 0.8 for x in xlim])
            self.ax.set_ylim([y * 0.8 for y in ylim])
            self.update_plot()

    def zoom_out(self):
        if self.ax:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self.ax.set_xlim([x * 1.2 for x in xlim])
            self.ax.set_ylim([y * 1.2 for y in ylim])
            self.update_plot()

    def reset_view(self):
        if self.ax:
            self.ax.set_xlim([-10, 10])
            self.ax.set_ylim([-10, 10])
            self.update_plot()

    def update_plot(self):
        if self.figure:
            self.figure.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = GraphingCalculator(root)
    root.mainloop()